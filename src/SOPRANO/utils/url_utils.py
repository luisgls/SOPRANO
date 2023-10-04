import pathlib
from pathlib import Path

import requests

from SOPRANO.utils.path_utils import Directories
from SOPRANO.utils.print_utils import task_output


class DownloadError(Exception):
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


def find_latest_release(partial_url: str):
    pass  # TODO


def build_ensembl_urls(species: str, reference: str):
    mixed_case_species = species[0].upper() + species[1:]
    base = (
        "https://ftp.ensembl.org/pub/release-{RELEASE}"
        + f"/fasta/{species}/dna"
    )
    partial = f"{base}/{mixed_case_species}.{reference}.dna"

    toplevel = f"{partial}.toplevel.fa.gz"
    primary_assembly = f"{partial}.primary_assembly.fa.gz"

    return {"toplevel": toplevel, "primary_assembly": primary_assembly}


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

    def download_toplevel(self, release=110):
        self._download(release, _toplevel=True)

    def download_primary_assembly(self, release=110):
        self._download(release, _toplevel=False)

    def decompress_toplevel(self, release=110):
        decompress(self.toplevel_fa_gz_path(release))

    def decompress_primary_assembly(self, release=110):
        decompress(self.primary_assembly_fa_gz_path(release))

    def compute_chrom_sizes(self, release=110):
        compute_chrom_sizes(self.primary_assembly_fa_path(release))


class EnsemblData(_GatherReferences):
    def __init__(self, species: str, reference: str):
        self.species = species
        self.reference = reference

        url_dict = build_ensembl_urls(species, reference)

        self.toplevel_url = url_dict["toplevel"]
        self.primary_assembly_url = url_dict["primary_assembly"]

    @classmethod
    def homo_sapiens_GRCh38(cls):
        return cls("homo_sapiens", "GRCh38")

    @classmethod
    def homo_sapiens_GRCh37(cls):
        return _GRCh37Data()


class _GRCh37Data(_GatherReferences):
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
