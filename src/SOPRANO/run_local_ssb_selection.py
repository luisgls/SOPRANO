import argparse
import pathlib
import subprocess

from SOPRANO import objects
from SOPRANO.misc_utils import Directories
from SOPRANO.pipeline_utils import run_pipeline


def check_path(cli_path: pathlib.Path | None, optional=False):
    if cli_path is None:
        if not optional:
            raise Exception("Input path is not optional and path is None!")
    elif not cli_path.exists():
        raise Exception(f"CLI input path does not exist: {cli_path}")


def check_genome(_namespace: argparse.Namespace) -> tuple:
    ref = _namespace.genome_ref
    release = _namespace.release

    available_refs = ("GRCh37", "GRCh38")

    if ref not in available_refs:
        raise ValueError(
            f"Reference {ref} not supported. "
            f"Permitted choices: {available_refs}"
        )

    return ref, release


def parse_args():
    parser = argparse.ArgumentParser(description="SOPRANO input arguments")

    parser.add_argument(
        "--input",
        "-i",
        dest="input_path",
        type=pathlib.Path,
        help="Prove the path to the input VEP annotated file.",
        required=True,
    )

    parser.add_argument(
        "--bed_file",
        "-b",
        dest="bed_path",
        type=pathlib.Path,
        help="Provide the path to the bed file with protein coordinates named "
        "by Transcript (ENSTXXXXXX 123 135)",
        required=True,
    )

    parser.add_argument(
        "--output",
        "-o",
        dest="cache_dir",
        type=pathlib.Path,
        help="Provide the path to the output directory in which dN/dS results "
        "will be cached.",
        required=True,
    )

    parser.add_argument(
        "--name",
        "-n",
        dest="analysis_name",
        type=str,
        help="Provide an identifying name for your results.",
        required=True,
    )
    analysis_params_group = parser.add_argument_group()

    analysis_params_group.add_argument(
        "--target_regions",
        "-t",
        dest="target_regions",
        type=pathlib.Path,
        help="Provide a bed file with regions to randomize.",
    )

    analysis_params_group.add_argument(
        "--use_ssb192",
        dest="use_ssb192",
        action="store_true",
        help="If flag is used ssb192 will be used, otherwise ssb7.",
    )

    analysis_params_group.add_argument(
        "--use_random",
        dest="use_random",
        action="store_true",
        help="If flag is used, calculate a dNdS calue for a random region "
        "similar to the target.",
    )

    # TODO: This is true by default in the original implementation
    analysis_params_group.add_argument(
        "--exclude_drivers",
        dest="exclude_drivers",
        action="store_true",
        help="If flag is used, driver geners will be excluded from the "
        "calculation.",
    )

    analysis_params_group.add_argument(
        "--seed",
        "-s",
        dest="seed",
        default=-1,
        type=int,
        help="Provide seed value for shuffle process in randomization. By "
        "default, seed value is < 0, for which, no seed value will be "
        "applied.",
    )

    transcript_args = parser.add_argument_group()

    transcript_args.add_argument(
        "--transcript",
        dest="transcript",
        help="Provide path to transcript file",
        default=objects.EnsemblTranscripts.transcript_length,
        type=pathlib.Path,
    )

    transcript_args.add_argument(
        "--protein_transcript",
        dest="protein_transcript",
        help="Provide path to protein transcript file",
        default=objects.EnsemblTranscripts.protein_transcript_length,
        type=pathlib.Path,
    )

    transcript_args.add_argument(
        "--fasta",
        dest="transcript_ids",
        help="Provide path to the ensembl transcript IDs fasta file",
        default=objects.EnsemblTranscripts.transcript_fasta,
        type=pathlib.Path,
    )

    genome_args = parser.add_argument_group()

    genome_args.add_argument(
        "--reference",
        "-r",
        dest="genome_ref",
        default="GRCh37",
        type=str,
        help="Reference genome file definition. By default, uses GRCh37. Pass "
        "instead GRCh38 if preferred. In order to download the reference "
        "fasta file, you can execute the command:\n"
        "GET_GENOMES -r GRCh37\n"
        "or\n"
        "GET_GENOMES -r GRCh38",
    )

    genome_args.add_argument(
        "--release",
        "-q",
        dest="release",
        type=str,
        help="Ensembl release number, e.g., 109, 110. Defaults to 110.",
        default="110",
    )

    args = parser.parse_args()

    check_path(args.input_path)
    check_path(args.bed_path)
    check_path(args.cache_dir)
    check_path(args.target_regions, optional=True)
    check_path(args.transcript)
    check_path(args.protein_transcript)
    check_path(args.transcript_ids)
    check_genome(args)

    return args


def title_output():
    print(
        """
███████  ██████  ██████  ██████   █████  ███    ██  ██████  
██      ██    ██ ██   ██ ██   ██ ██   ██ ████   ██ ██    ██ 
███████ ██    ██ ██████  ██████  ███████ ██ ██  ██ ██    ██ 
     ██ ██    ██ ██      ██   ██ ██   ██ ██  ██ ██ ██    ██ 
███████  ██████  ██      ██   ██ ██   ██ ██   ████  ██████  
"""
    )


def line_output(n=60):
    print("-" * n)


def startup_output(**kwargs):
    title_output()
    line_output()
    print("Selection On PRotein ANnotated regiOns")
    line_output()

    if kwargs:
        # Parameters used in pipeline
        for k, v in kwargs.items():
            print("-> {0:.<30}".format(k) + f"{v}")

        line_output()


def main(_namespace=None):
    cli_args = parse_args() if _namespace is None else _namespace
    startup_output(**cli_args.__dict__)
    params = objects.Parameters.from_namespace(cli_args)
    run_pipeline(params)


def parse_genome_args():
    parser = argparse.ArgumentParser(description="Genome reference")

    parser.add_argument(
        "--reference",
        "-r",
        dest="genome_ref",
        type=str,
        help="GRCh37 or GRCh38",
        required=True,
    )

    parser.add_argument(
        "--release",
        dest="release",
        type=str,
        help="Ensemblv release number, e.g., 109, 110. Defaults to 110.",
        default="110",
    )

    ref_release = check_genome(parser.parse_args())

    return ref_release


def download_genome():
    ref, release = parse_genome_args()
    startup_output()
    downloader_path = Directories.installers("download_reference.sh")
    data_dir = Directories.homo_sapien_genomes(f"{release}_{ref}")

    if not data_dir.exists():
        data_dir.mkdir(parents=True)

    assert downloader_path.exists(), downloader_path
    assert data_dir.exists(), data_dir

    subprocess.run(
        ["bash", downloader_path.as_posix(), ref, release, data_dir.as_posix()]
    )


def local_st_app():
    """
    Runs streamlit app interface for SOPRANO
    """
    app_path = Directories.src("st_app.py")
    subprocess.run(["streamlit", "run", app_path])


if __name__ == "__main__":
    main()
