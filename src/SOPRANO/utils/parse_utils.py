import argparse
import pathlib

from SOPRANO.core import objects
from SOPRANO.utils.path_utils import Directories, check_cli_path


def _add_core_genome_args(parser: argparse.ArgumentParser):
    parser.add_argument(
        "--species",
        "-s",
        dest="species",
        type=str,
        help="Ensembl species Latin name. " "Default: homo_sapiens.",
        default="homo_sapiens",
    )

    parser.add_argument(
        "--assembly",
        "-a",
        dest="assembly",
        type=str,
        help="Ensembl genome assembly ID. Default: GRCh38.",
        default="GRCh38",
    )

    parser.add_argument(
        "--release",
        "-r",
        dest="release",
        type=int,
        help="Ensembl release number. Default: 110.",
        default=110,
    )

    return parser


def fix_species_arg(species: str) -> str:
    return species.replace(" ", "_").lower()


def fix_ns_species_arg(_namespace: argparse.Namespace) -> argparse.Namespace:
    _namespace.species = fix_species_arg(_namespace.species)
    return _namespace


def parse_genome_args(argv=None):
    parser = argparse.ArgumentParser(description="Genome reference")
    parser = _add_core_genome_args(parser)

    parser.add_argument(
        "--primary_assembly",
        "-p",
        dest="primary_assembly",
        action="store_true",
    )

    parser.add_argument(
        "--download_only", "-d", dest="download_only", action="store_true"
    )
    return fix_ns_species_arg(parser.parse_args(argv))


def parse_args(argv=None):
    parser = argparse.ArgumentParser(description="SOPRANO input arguments")

    parser.add_argument(
        "--input",
        "-i",
        dest="input_path",
        type=pathlib.Path,
        help="Prove the path to the input annotated file.",
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
        dest="random_regions",
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
        "-x",
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
        help=f"Provide path to transcript file. Default: "
        f"{objects.TranscriptPaths.defaults().transcript_length}",
        default=objects.TranscriptPaths.defaults().transcript_length,
        type=pathlib.Path,
    )

    transcript_args.add_argument(
        "--protein_transcript",
        "-p",
        dest="protein_transcript",
        help=f"Provide path to protein transcript file. Default: "
        f"{objects.TranscriptPaths.defaults().protein_transcript_length}",
        default=objects.TranscriptPaths.defaults().protein_transcript_length,
        type=pathlib.Path,
    )

    transcript_args.add_argument(
        "--fasta",
        "-f",
        dest="transcript_ids",
        help=f"Provide path to the ensembl transcript IDs fasta file. "
        f"Default: {objects.TranscriptPaths.defaults().transcript_fasta}",
        default=objects.TranscriptPaths.defaults().transcript_fasta,
        type=pathlib.Path,
    )

    _add_core_genome_args(parser)

    args = parser.parse_args(argv)
    args = fix_ns_species_arg(args)

    check_cli_path(args.input_path)
    check_cli_path(args.bed_path)
    check_cli_path(args.cache_dir)
    check_cli_path(args.random_regions, optional=True)
    check_cli_path(args.transcript)
    check_cli_path(args.protein_transcript)
    check_cli_path(args.transcript_ids)

    return args


def parse_hla(argv=None):
    parser = argparse.ArgumentParser("Parse HLA parameters")
    parser.add_argument(
        "--alleles",
        "-a",
        dest="hla_values",
        nargs="+",
        type=str,
        required=True,
        help="Space seperated HLA alleles.",
    )
    parser.add_argument(
        "--output",
        "-o",
        dest="output_id",
        type=str,
        required=True,
        help="Identifying name for the output immunopeptidome file.",
    )
    parser.add_argument(
        "--cache",
        "-c",
        dest="cache_dir",
        type=pathlib.Path,
        default=Directories.app_immunopeptidomes(),
        help=f"Cache location for immunopeptidome file. "
        f"Default: {Directories.app_immunopeptidomes()}",
    )
    parser.add_argument(
        "--restrict",
        "-r",
        dest="restricted_transcript_ids",
        nargs="*",
        type=str,
        required=False,
        help="Space seperated Ensembl transcript IDs. Default: []",
        default=[],
    )
    parser.add_argument(
        "--excluded",
        "-e",
        dest="excluded_transcript_ids",
        nargs="*",
        type=str,
        required=False,
        help="Space seperated Ensembl transcript IDs. Default: []",
        default=[],
    )

    return parser.parse_args(argv)


def parse_vcf_sources(argv=None):
    parser = argparse.ArgumentParser(description="VCF file annotations")

    parser.add_argument(
        "--source_dir",
        "-s",
        dest="sources_dir",
        type=pathlib.Path,
        help="Provide path to directory containing VCF files.",
    )
    parser.add_argument(
        "--output",
        "-o",
        dest="output_name",
        type=str,
        help="Provide name for output file. No extension required.",
    )
    parser.add_argument(
        "--cache_dir",
        "-d",
        dest="cache_dir",
        type=pathlib.Path,
        default=Directories.app_annotated_inputs(),
        help=f"Provide path to directory that will contain output. Default: "
        f"{Directories.app_annotated_inputs()}",
    )

    return parser.parse_args(argv)
