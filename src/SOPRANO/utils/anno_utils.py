from pathlib import Path

from SOPRANO.utils.path_utils import Directories
from SOPRANO.utils.sh_utils import pipe


def find_vcf_files(vcf_sources_dir: Path):
    vcf_exts = {".vcf", ".VCF"}
    gz_exts = {".gz", ".GZ", ".Gz"}

    exts = {x + y for x in vcf_exts for y in gz_exts}

    detected = []

    for ext in exts:
        for vcf_file_path in vcf_sources_dir.glob(f"*{ext}"):
            detected.append(vcf_file_path)

    n_detected = len(detected)

    if n_detected > 0:
        print(
            f"Detected {n_detected} vcf files in {vcf_sources_dir.as_posix()}"
        )
    else:
        raise FileNotFoundError(
            f"No VCF files detected in {vcf_sources_dir.as_posix()}"
        )

    return detected


def annotate_source(
    vcf_sources_dir: Path,
    output_name: str | None = None,
    cache_directory: Path = Directories.app_annotated_inputs(),
):
    rscript_path = Directories.r_scripts("parse_vcf.R")

    if output_name is None:
        output_path = cache_directory.joinpath(vcf_sources_dir.name)
    else:
        output_path = cache_directory.joinpath(output_name)

    output_path = output_path.with_suffix(".anno")

    find_vcf_files(vcf_sources_dir)

    pipe(
        [
            "Rscript",
            rscript_path.as_posix(),
            "-o",
            output_path.as_posix(),
            "-d",
            vcf_sources_dir.as_posix(),
            "-t",
            Directories.data().as_posix(),
        ]
    )


if __name__ == "__main__":
    x = Path("/mnt/c/Users/kmarzouk/software/misc/vcf_parser/data")
    annotate_source(x)
