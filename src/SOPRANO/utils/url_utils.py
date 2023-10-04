import pathlib
from pathlib import Path

import requests

from SOPRANO.utils.path_utils import Directories
from SOPRANO.utils.print_utils import task_output


class DownloadError(Exception):
    pass


def filename_from_url(url: str):
    return url.split("/")[-1]


def download_from_url(url: str, target_dir: Path = Path.cwd()):
    task_output(f"Attempting request from: {url}")
    response = requests.get(url)

    if response.status_code != 200:
        raise DownloadError(f"Download status code = {response.status_code}")

    if "content-disposition" in response.headers:
        content_disposition = response.headers["content-disposition"]
        filename = content_disposition.split("filename=")[1]
    else:
        filename = filename_from_url(url)

    download_path = target_dir.joinpath(filename)

    task_output(f"Writing content to {download_path.as_posix()}")
    with open(download_path, mode="wb") as file:
        file.write(response.content)
    task_output("Process complete.")


def decompress(gz_path: pathlib.Path):
    pass  # TODO


def compute_chrom_count(path: pathlib.Path):
    pass  # TODO


def find_latest_release(partial_url: str):
    pass  # TODO


class _GetGenome:
    # Urls
    toplevel_url: str
    primary_assembly_url: str

    # Params
    species: str
    reference: str

    # Status
    toplevel_done: bool
    primary_assembly_done: bool
    sizes_done: bool

    def _download(self, release: int, _toplevel):
        species_dir = Directories.data(self.species)
        release_dir = species_dir.joinpath(self.species).joinpath(
            f"{release}_{self.reference}"
        )

        url = self.toplevel_url if _toplevel else self.primary_assembly_url

        download_from_url(url.format(release=release), target_dir=release_dir)

    def download_toplevel(self, release=110):
        self._download(release, _toplevel=True)

    def download_primary_assembly(self, release=110):
        self._download(release, _toplevel=False)
