import gzip
import pathlib
import pickle as pk
import shutil
from pathlib import Path

import requests
from clint.textui import progress

from SOPRANO.utils.path_utils import Directories
from SOPRANO.utils.print_utils import task_output
from SOPRANO.utils.sh_utils import pipe

_SOPRANO_ENSEMBL_RELEASES = Directories.data("ensembl.releases")


class DownloadError(Exception):
    pass


class BadResponse(Exception):
    pass


def filename_from_url(url: str):
    return url.split("/")[-1]


def download_from_url(url: str, target_path: Path | None = None):
    task_output(f"Attempting request from: {url}")

    response = requests.get(url, stream=True)

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

    with open(target_path, "wb") as f:
        # fmt: off
        total_length = int(response.headers.get("content-length"))  # type: ignore[arg-type]
        # fmt: on
        for chunk in progress.bar(
            response.iter_content(chunk_size=1024),
            expected_size=(total_length / 1024) + 1,
        ):
            if chunk:
                f.write(chunk)
                f.flush()


def decompress(gz_path: pathlib.Path):
    if not gz_path.exists():
        raise FileNotFoundError(gz_path)

    output_path = gz_path.with_suffix("")

    task_output(f"Decompressing {gz_path}")
    with gzip.open(gz_path, "rb") as f_in:
        with open(output_path, "wb") as f_out:
            shutil.copyfileobj(f_in, f_out)


def compute_fasta_index(fa_path: pathlib.Path):
    if not fa_path.exists():
        raise FileNotFoundError(fa_path)

    index_path = fa_path.with_suffix(".fa.fai")

    task_output(f"Computing fasta index file from {fa_path}")
    pipe(
        ["samtools", "faidx", fa_path.as_posix(), "-o", index_path.as_posix()]
    )


def compute_chrom_sizes(fai_path: pathlib.Path):
    if not fai_path.exists():
        raise FileNotFoundError(fai_path)

    if not fai_path.name.endswith(".fa.fai"):
        raise ValueError(f"Bad extension: {fai_path}")

    chrom_path = fai_path.with_suffix("").with_suffix(".chrom")

    task_output(f"Computing chromosome sizes from {fai_path}")
    pipe(["cut", "-f1,2", fai_path.as_posix()], output_path=chrom_path)


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


def check_ensembl_file_url(file_url: str, release: int):
    file_url = file_url.format(RELEASE=release)
    response = requests.head(file_url)
    if response.status_code != 200:
        raise BadResponse(
            f"Response for {file_url} = {response.status_code}\nIt is "
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
            check_ensembl_file_url(toplevel_url, release)
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
            check_ensembl_file_url(toplevel_url, release)
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
