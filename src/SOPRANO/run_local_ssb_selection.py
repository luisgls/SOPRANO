import argparse
import pathlib
from datetime import datetime

from SOPRANO import objects, prepare_coordinates


def check_path(cli_path: pathlib.Path | None, optional=False):
    if cli_path is None:
        if not optional:
            raise Exception("Input path is not optional and path is None!")
    elif not cli_path.exists():
        raise Exception(f"CLI input path does not exist: {cli_path}")


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

    analysis_params_group.add_argument(
        "--exclude_drivers",
        dest="exclude_drivers",
        action="store_true",
        help="If flag is used, driver geners will be excluded from the "
        "calculation.",
    )

    analysis_params_group.add_argument(
        "--seed",
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

    args = parser.parse_args()

    check_path(args.input_path)
    check_path(args.bed_path)
    check_path(args.cache_dir)
    check_path(args.target_regions, optional=True)
    check_path(args.transcript)
    check_path(args.protein_transcript)

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


def time_output():
    now = datetime.now()
    return now.strftime("%d/%m/%Y %H:%M:%S")


def task_output(msg):
    print(f"[{time_output()}] {msg}")


def line_output(n=60):
    print("-" * n)


def startup_output(**kwargs):
    title_output()
    line_output()
    print("Selection On PRotein ANnotated regiOns")
    line_output()

    # Parameters used in pipeline
    for k, v in kwargs.items():
        print("-> {0:.<30}".format(k) + f"{v}")

    line_output()


def main(_namespace=None):
    cli_args = parse_args() if _namespace is None else _namespace
    startup_output(**cli_args.__dict__)

    params = objects.Parameters.from_namespace(cli_args)

    task_output("Filtering transcripts")
    prepare_coordinates.FilterTranscripts.apply(params)

    if params.use_target_regions:
        task_output(
            f"Randomizing transcripts subject to bed file: "
            f"{params.target_regions_path.as_posix()}"
        )
        randomization_method = prepare_coordinates.RandomizeWithRegions
    elif params.use_random:
        task_output("Randomizing transcripts")
        randomization_method = prepare_coordinates.RandomizeWithoutRegions
    else:
        task_output(
            f"No randomization selected: sorting input bed file: "
            f"{params.bed_path.as_posix()}"
        )
        randomization_method = prepare_coordinates.NonRandom

    randomization_method.prep_epitopes_ori2(params)

    if params.exclude_drivers:
        task_output("Excluding positively selected genes")
        drivers_method = prepare_coordinates.GeneExclusions
    else:
        task_output("Retaining positively selected genes")
        drivers_method = prepare_coordinates.GeneExclusionsDisabled

    drivers_method.apply(params)


if __name__ == "__main__":
    main()
