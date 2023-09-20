"""

vep -i homo_sapiens_GRCh38.vcf --cache --force_overwrite --dir_cache ../data/
    --assembly GRCh38

seems to work from the examples dir ... but not in format we want

"""

import pathlib

from SOPRANO.sh_utils import subprocess_pipes

_DIR_CACHE = pathlib.Path(__file__).parent.joinpath("data")
_EXAMPLES = pathlib.Path(__file__).parent.joinpath("examples")
_FASTA = _DIR_CACHE.joinpath("ensemble_transcriptID.translated.fasta")


def annotate_homo_sapiens(
    input_file: pathlib.Path,
    assembly="GRCh38",
    output_file: pathlib.Path | None = None,
    overwrite=False,
):
    if output_file is None:
        output_file = input_file.with_suffix(".anno")

    subprocess_pipes.pipe(
        [
            "vep",
            "-i",
            input_file.as_posix(),
            "-o",
            output_file.as_posix(),
            "--cache",
            "--dir_cache",
            _DIR_CACHE.as_posix(),
            "--assembly",
            assembly,
            "--all_refseq",
            "--pick",
            "--symbol",
            "--no_stats",
            "--force_overwrite",  # TODO: Fix optional
            "--fasta",
            _FASTA.as_posix(),
        ]
    )


input_vcf = _EXAMPLES.joinpath("homo_sapiens_GRCh38.vcf")
output = _EXAMPLES.joinpath("with_fasta.anno")
annotate_homo_sapiens(input_vcf, output_file=output)
