import argparse
import pathlib
from datetime import datetime

from SOPRANO import objects


def _get_datetime_str():
    now = datetime.now()
    return now.strftime("%d/%m/%Y %H:%M:%S")


def task_output(msg):
    print(f"[{_get_datetime_str()}] {msg}")


def parse_args():
    parser = argparse.ArgumentParser(description="SOPRANO input arguments")

    parser.add_argument(
        "--input",
        "-i",
        dest="input",
        type=pathlib.Path,
        help="Prove the path to the input VEP annotated file.",
        required=True,
    )

    parser.add_argument(
        "--bed_file",
        "-b",
        dest="bed_file",
        type=pathlib.Path,
        help="Provide the path to the bed file with protein coordinates named "
        "by Transcript (ENSTXXXXXX 123 135)",
        required=True,
    )

    parser.add_argument(
        "--output",
        "-o",
        dest="output",
        type=pathlib.Path,
        help="Provide the path to the output directory in which dN/dS results "
        "will be cached.",
        required=True,
    )

    parser.add_argument(
        "--name",
        "-n",
        dest="name",
        type=str,
        help="Provide an identifying name for your results.",
        required=True,
    )
    analysis_params_group = parser.add_argument_group()

    analysis_params_group.add_argument(
        "-t",
        dest="bed_regions",
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
        "--exclude_drivers",
        action="store_true",
        help="If flag is used, driver geners will be excluded from the "
        "calculation.",
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

    args = parser.parse_args()

    _validate_input_path(args.input)
    _validate_input_path(args.bed_file)
    _validate_input_path(args.output)
    _validate_input_path(args.bed_regions, optional=True)
    _validate_input_path(args.transcript)
    _validate_input_path(args.protein_transcript)

    return args


def _validate_input_path(cli_path: pathlib.Path | None, optional=False):
    if cli_path is None:
        if not optional:
            raise Exception("Input path is not optional and path is None!")
    elif not cli_path.exists():
        raise Exception(f"CLI input path does not exist: {cli_path}")


_title = """
███████  ██████  ██████  ██████   █████  ███    ██  ██████  
██      ██    ██ ██   ██ ██   ██ ██   ██ ████   ██ ██    ██ 
███████ ██    ██ ██████  ██████  ███████ ██ ██  ██ ██    ██ 
     ██ ██    ██ ██      ██   ██ ██   ██ ██  ██ ██ ██    ██ 
███████  ██████  ██      ██   ██ ██   ██ ██   ████  ██████  
"""
_borders = "-" * 60


def _startup_message(**kwargs):
    print(_title)
    print(_borders)
    print("Selection On PRotein ANnotated regiOns")
    print(_borders)

    for k, v in kwargs.items():
        print("-> {0:.<30}".format(k) + f"{v}")

    print(_borders)


def main(_namespace=None):
    cli_args = parse_args() if _namespace is None else _namespace
    _startup_message(**cli_args.__dict__)

    task_output("Filtering transcripts")

    if cli_args.bed_regions is None:
        bed_regions = None
    else:
        bed_regions = pathlib.Path(cli_args.bed_regions)

    transcripts = objects.TranscriptPaths(
        pathlib.Path(cli_args.transcript),
        pathlib.Path(cli_args.protein_transcript),
    )

    params = objects.Parameters(
        analysis_name=cli_args.name,
        bed_path=cli_args.bed_file,
        tmpdir=pathlib.Path(cli_args.output).joinpath("tmp"),
        target_regions_path=bed_regions,
        transcripts=transcripts,
    )


if __name__ == "__main__":
    main()
