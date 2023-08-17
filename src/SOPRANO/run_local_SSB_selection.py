import argparse


def main(*args, **kwargs):
    pass


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="SOPRANO input arguments")

    parser.add_argument(
        "--bed_file",
        "-b",
        dest="bed_file",
        type=str,
        help="Provide the path to the bed file with protein coordinates named "
             "by Transcript (ENSTXXXXXX 123 135)",
        required=True
    )

    parser.add_argument(
        "--output",
        "-o",
        dest="output",
        type=str,
        help="Provide the path to the output directory in which dN/dS results "
             "will be cached.",
        required=True
    )

    parser.add_argument(
        "--name",
        "-n",
        dest="name",
        type=str,
        help="Provide an identifying name for your results.",
        required=True
    )

    args = parser.parse_args()
