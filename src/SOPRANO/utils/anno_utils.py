from pathlib import Path

from SOPRANO.utils.path_utils import Directories
from SOPRANO.utils.sh_utils import pipe


class NotGZ(Exception):
    pass


class NotVCF(Exception):
    pass


def find_vcf_files(vcf_source: Path):
    vcf_exts = {".vcf", ".VCF"}
    gz_exts = {".gz", ".GZ", ".Gz"}

    exts = {x + y for x in vcf_exts for y in gz_exts}

    if vcf_source.is_file():
        if vcf_source.suffix not in gz_exts:
            raise NotGZ(
                f"Input source should be gzip compressed: " f"{vcf_source}"
            )

        if (vcf_source.with_suffix("")).suffix not in vcf_exts:
            raise NotVCF(
                f"Input sourse should contain vcf extension: " f"{vcf_source}"
            )

        return [vcf_source]
    elif not vcf_source.is_dir():
        raise ValueError(
            "Sources input must be a file of vcf.gz format, or a directory"
            f"containing vcf.gz files. {vcf_source} is invalid."
        )

    detected = []

    for ext in exts:
        for vcf_file_path in vcf_source.glob(f"*{ext}"):
            detected.append(vcf_file_path)

    n_detected = len(detected)

    if n_detected > 0:
        print(f"Detected {n_detected} vcf files in {vcf_source.as_posix()}")
    else:
        raise FileNotFoundError(
            f"No VCF files detected in {vcf_source.as_posix()}"
        )

    return detected


def annotate_source(
    vcf_sources_dir: Path,
    assembly: str,
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
            "-a",
            assembly,
        ]
    )
