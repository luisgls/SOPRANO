import argparse
import pathlib


def main(*args, **kwargs):
    pass


def _validate_input_path(cli_path: pathlib.PosixPath, optional=False):
    pass


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="SOPRANO input arguments")

    parser.add_argument(
        "--input",
        "-i",
        dest="input",
        type=pathlib.PosixPath,
        help="Prove the path to the input VEP annotated file.",
        required=True,
    )

    parser.add_argument(
        "--bed_file",
        "-b",
        dest="bed_file",
        type=pathlib.PosixPath,
        help="Provide the path to the bed file with protein coordinates named "
        "by Transcript (ENSTXXXXXX 123 135)",
        required=True,
    )

    parser.add_argument(
        "--output",
        "-o",
        dest="output",
        type=pathlib.PosixPath,
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
        type=pathlib.PosixPath,
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

    args = parser.parse_args()

    _validate_input_path(args.input)
    _validate_input_path(args.bed_file)
    _validate_input_path(args.output)
    _validate_input_path(args.bed_regions, optional=True)
