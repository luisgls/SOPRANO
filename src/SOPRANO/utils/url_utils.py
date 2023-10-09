import pathlib
import pickle as pk
from pathlib import Path

import requests

from SOPRANO.utils.path_utils import Directories
from SOPRANO.utils.print_utils import task_output

_SOPRANO_ENSEMBL_RELEASES = Directories.data("ensembl.releases")


class DownloadError(Exception):
    pass


class BadResponse(Exception):
    pass


def filename_from_url(url: str):
    return url.split("/")[-1]


def download_from_url(url: str, target_path: Path | None = None):
    task_output(f"Attempting request from: {url}")
    response = requests.get(url)

    if response.status_code != 200:
        raise DownloadError(f"Download status code = {response.status_code}")

    if "content-disposition" in response.headers:
        content_disposition = response.headers["content-disposition"]
        filename = content_disposition.split("filename=")[1]
    else:
        filename = filename_from_url(url)

    if target_path is None:
        target_path = Path.cwd().joinpath(filename)

    task_output(f"Writing content to {target_path.as_posix()}")
    with open(target_path, mode="wb") as file:
        file.write(response.content)
    task_output("Process complete.")


def decompress(gz_path: pathlib.Path):
    pass  # TODO


def compute_chrom_sizes(path: pathlib.Path):
    pass  # TODO


def get_dna_url(species: str):
    return (
        "https://ftp.ensembl.org/pub/release-{RELEASE}"
        + f"/fasta/{species}/dna"
    )


def check_ensembl_dna_url(dna_url: str, release: int):
    dna_url = dna_url.format(RELEASE=release)
    response = requests.get(dna_url)
    if response.status_code != 200:
        raise BadResponse(
            f"Response for {dna_url} = {response.status_code}\n"
            f"It is likely that "
            f"there is no ensembl release {release} for the target genome "
            f"reference."
        )


def check_ensembl_toplevel_url(toplevel_url: str, release: int):
    toplevel_url = toplevel_url.format(RELEASE=release)
    response = requests.head(toplevel_url)
    if response.status_code != 200:
        raise BadResponse(
            f"Response for {toplevel_url} = {response.status_code}\nIt is "
            f"likely that there is no ensembl release {release} for the "
            f"target genome reference."
        )


def get_cached_release_dict():
    try:
        with open(_SOPRANO_ENSEMBL_RELEASES, "rb") as f:
            release_dict = pk.load(f)
    except FileNotFoundError:
        release_dict = {}

    return release_dict


def get_cached_release_value(dna_url: str):
    release_dict = get_cached_release_dict()
    return release_dict.get(dna_url, [-1, -1])


def _cache_release_data(toplevel_url: str, release: int, _earliest: bool):
    cache_index = 0 if _earliest else 1
    cached_release_dict = get_cached_release_dict()
    releases_list = cached_release_dict.get(toplevel_url, [-1, -1])

    if releases_list[cache_index] < release:
        releases_list[cache_index] = release

    cached_release_dict[toplevel_url] = releases_list

    _SOPRANO_ENSEMBL_RELEASES.unlink(missing_ok=True)
    with open(_SOPRANO_ENSEMBL_RELEASES, "wb") as f:
        pk.dump(cached_release_dict, f)


def cache_earliest_release(toplevel_url: str, release: int):
    _cache_release_data(toplevel_url, release, True)


def cache_latest_release(toplevel_url: str, release: int):
    _cache_release_data(toplevel_url, release, False)


def find_earliest_release(toplevel_url: str, _max_release=200):
    cached_earliest = get_cached_release_value(toplevel_url)[0]

    if cached_earliest > -1:
        return cached_earliest

    release = 0

    while release < _max_release:
        try:
            check_ensembl_toplevel_url(toplevel_url, release)
            break
        except BadResponse:
            release += 1

    if release == _max_release:
        raise BadResponse(f"Unable to find release content for {toplevel_url}")

    cache_earliest_release(toplevel_url, release)

    return release


def find_latest_release(toplevel_url: str, _max_release=200):
    release = max(get_cached_release_value(toplevel_url))

    if release == -1:
        release = find_earliest_release(toplevel_url, _max_release)

    while release < _max_release:
        try:
            check_ensembl_toplevel_url(toplevel_url, release)
            release += 1
        except BadResponse:
            release -= 1
            break

    cache_latest_release(toplevel_url, release)

    return release


def build_ensembl_urls(species: str, reference: str):
    mixed_case_species = species[0].upper() + species[1:]
    dna_url = get_dna_url(species)
    _url = f"{dna_url}/{mixed_case_species}.{reference}.dna"

    toplevel_url = f"{_url}.toplevel.fa.gz"
    primary_assembly_url = f"{_url}.primary_assembly.fa.gz"

    return {"toplevel": toplevel_url, "primary_assembly": primary_assembly_url}


class _GatherReferences:
    # Urls
    toplevel_url: str
    primary_assembly_url: str

    # Params
    species: str
    reference: str

    # Status
    toplevel_done: bool = False
    primary_assembly_done: bool = False
    sizes_done: bool = False

    def _dest_directory(self, release: int):
        return Directories.data(self.species) / f"{release}_{self.reference}"

    def _dest_fa_gz(self, release: int, _toplevel: bool):
        return self._dest_directory(release) / filename_from_url(
            self.toplevel_url if _toplevel else self.primary_assembly_url
        ).format(RELEASE=release)

    def _dest_fa(self, release: int, _toplevel: bool):
        return self._dest_fa_gz(release, _toplevel).with_suffix("")

    def dest_chrom(self, release: int, _toplevel: bool):
        return self._dest_fa(release, _toplevel).with_suffix(".chrom")

    def toplevel_fa_gz_path(self, release: int):
        return self._dest_fa_gz(release, _toplevel=True)

    def toplevel_fa_path(self, release: int):
        return self._dest_fa(release, _toplevel=True)

    def toplevel_chrom_path(self, release: int):
        return self.toplevel_fa_path(release).with_suffix(".chrom")

    def primary_assembly_fa_gz_path(self, release: int):
        return self._dest_fa_gz(release, _toplevel=False)

    def primary_assembly_fa_path(self, release: int):
        return self._dest_fa(release, _toplevel=False)

    def _download(self, release: int, _toplevel):
        if _toplevel:
            source_url = self.toplevel_url
            dest_path = self._dest_fa_gz(release, _toplevel)
        else:
            source_url = self.primary_assembly_url
            dest_path = self._dest_fa_gz(release, _toplevel)

        download_from_url(
            source_url,
            target_path=dest_path,
        )

    def _check_release_ok(self, release):
        min_release = find_earliest_release(self.toplevel_url)
        max_release = find_latest_release(self.toplevel_url)

        if not (min_release <= release <= max_release):
            raise ValueError(release)

    def download_toplevel(self, release=110):
        self._check_release_ok(release)
        self._download(release, _toplevel=True)

    def download_primary_assembly(self, release=110):
        self._check_release_ok(release)
        self._download(release, _toplevel=False)

    def decompress_toplevel(self, release=110):
        decompress(self.toplevel_fa_gz_path(release))

    def decompress_primary_assembly(self, release=110):
        decompress(self.primary_assembly_fa_gz_path(release))

    def compute_chrom_sizes(self, release=110):
        compute_chrom_sizes(self.primary_assembly_fa_path(release))


class EnsemblData(_GatherReferences):
    def __init__(self, species: str, reference: str, _init_urls=True):
        self.species = species
        self.reference = reference

        if _init_urls:
            url_dict = build_ensembl_urls(species, reference)
            self.toplevel_url = url_dict["toplevel"]
            self.primary_assembly_url = url_dict["primary_assembly"]

    @classmethod
    def homo_sapiens_GRCh38(cls):
        return cls("homo_sapiens", "GRCh38")

    @classmethod
    def homo_sapiens_GRCh37(cls):
        # GRCh37 is actually has a deviant url structure, so manually set here
        toplevel_url = (
            "https://ftp.ensembl.org/pub/grch37/release-{RELEASE}/"
            "fasta/homo_sapiens/dna/"
            "Homo_sapiens.GRCh37.dna.toplevel.fa.gz"
        )

        primary_assembly_url = (
            "https://ftp.ensembl.org/pub/grch37/release-{RELEASE}/"
            "fasta/homo_sapiens/dna/"
            "Homo_sapiens.GRCh37.dna.primary_assembly.fa.gz"
        )

        species = "homo_sapiens"
        reference = "GRCh37"

        obj = cls(species, reference, _init_urls=False)
        obj.toplevel_url = toplevel_url
        obj.primary_assembly_url = primary_assembly_url
        return obj
