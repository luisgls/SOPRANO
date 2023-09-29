import argparse
import pathlib

from SOPRANO.core import objects
from SOPRANO.utils.misc_utils import check_cli_path


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
        "--random_regions",
        "-m",
        dest="random_regions",  # TODO: Update to random_regions
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

    analysis_params_group.add_argument(
        "--keep_drivers",
        dest="keep_drivers",
        action="store_true",
        help="If flag is used, driver genes in src/data/genes2exclude.txt will"
        " be retained in the calculation.",
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
        "-t",
        dest="transcript",
        help="Provide path to transcript file",
        default=objects.EnsemblTranscripts.transcript_length,
        type=pathlib.Path,
    )

    transcript_args.add_argument(
        "--protein_transcript",
        "-p",
        dest="protein_transcript",
        help="Provide path to protein transcript file",
        default=objects.EnsemblTranscripts.protein_transcript_length,
        type=pathlib.Path,
    )

    transcript_args.add_argument(
        "--fasta",
        "-f",
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

    check_cli_path(args.input_path)
    check_cli_path(args.bed_path)
    check_cli_path(args.cache_dir)
    check_cli_path(args.random_regions, optional=True)
    check_cli_path(args.transcript)
    check_cli_path(args.protein_transcript)
    check_cli_path(args.transcript_ids)
    check_genome(args)

    return args
