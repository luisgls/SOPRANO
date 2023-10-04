from pathlib import Path

import requests

from SOPRANO.utils.print_utils import task_output


class DownloadError(Exception):
    pass


def download_file(url: str, target_dir: Path = Path.cwd()):
    task_output(f"Attempting request from: {url}")
    response = requests.get(url)

    if response.status_code != 200:
        raise DownloadError(f"Download status code = {response.status_code}")

    if "content-disposition" in response.headers:
        content_disposition = response.headers["content-disposition"]
        filename = content_disposition.split("filename=")[1]
    else:
        filename = url.split("/")[-1]

    download_path = target_dir.joinpath(filename)

    task_output(f"Writing content to {download_path.as_posix()}")
    with open(download_path, mode="wb") as file:
        file.write(response.content)
    task_output("Process complete.")
