import subprocess

from SOPRANO.core import objects
from SOPRANO.pipeline import run_pipeline
from SOPRANO.utils.misc_utils import Directories
from SOPRANO.utils.parse_utils import parse_args, parse_genome_args
from SOPRANO.utils.print_utils import startup_output
from SOPRANO.utils.vep_utils import (
    _get_src_dst_link_pairs,
    _link_src_dst_pairs,
    _link_vep_cache_parser,
)


def run_cli(_namespace=None):
    cli_args = parse_args() if _namespace is None else _namespace
    startup_output(**cli_args.__dict__)
    params = objects.Parameters.from_namespace(cli_args)
    run_pipeline(params)


def run_app():
    app_path = Directories.src("app.py")
    subprocess.run(["streamlit", "run", app_path])


def link_vep_cache():
    src_cache = _link_vep_cache_parser()
    src_dst_links = _get_src_dst_link_pairs(src_cache)
    _link_src_dst_pairs(src_dst_links)


def download_genome():
    ref, release = parse_genome_args()
    startup_output()
    downloader_path = Directories.installers("download_homo_sapiens.sh")
    data_dir = Directories.homo_sapien_genomes(f"{release}_{ref}")

    if not data_dir.exists():
        data_dir.mkdir(parents=True)

    assert downloader_path.exists(), downloader_path
    assert data_dir.exists(), data_dir

    subprocess.run(
        ["bash", downloader_path.as_posix(), ref, release, data_dir.as_posix()]
    )


if __name__ == "__main__":
    run_cli()
