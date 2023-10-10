import subprocess

from SOPRANO import hla2ip
from SOPRANO.core import objects
from SOPRANO.pipeline import run_pipeline
from SOPRANO.utils.parse_utils import parse_args, parse_genome_args, parse_hla
from SOPRANO.utils.path_utils import Directories
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
    subprocess.run(["streamlit", "run", app_path.as_posix()])


def link_vep_cache():
    src_cache = _link_vep_cache_parser()
    src_dst_links = _get_src_dst_link_pairs(src_cache)
    _link_src_dst_pairs(src_dst_links)


def download_genome():
    args = parse_genome_args()
    startup_output(**args.__dict__)
    print("Starting download...")

    if args.assembly == "GRCh37":
        assert args.species == "homo_sapiens"
        ensembl_data = objects.EnsemblData.homo_sapiens_GRCh37()
    else:
        ensembl_data = objects.EnsemblData(args.species, args.assembly)

    if args.primary_assembly:
        ensembl_data.download_primary_assembly(args.release)
        if not args.download_only:
            ensembl_data.compute_all_primary_assembly(args.release)
    else:
        ensembl_data.download_toplevel(args.release)
        if not args.download_only:
            ensembl_data.compute_all_toplevel(args.release)


def hla2pip():
    args = parse_hla()
    hla2ip.immunopeptidome_from_hla(
        *args.hla_values, output_name=args.output_id, cache_loc=args.cache_dir
    )


if __name__ == "__main__":
    run_cli()
