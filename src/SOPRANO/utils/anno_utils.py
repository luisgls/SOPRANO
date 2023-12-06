from pathlib import Path

from SOPRANO.utils.path_utils import Directories
from SOPRANO.utils.sh_utils import pipe


class NotVCF(Exception):
    pass


class NoVCFs(Exception):
    pass


def find_vcf_files(vcf_source: Path):
    vcf_exts = {".vcf", ".VCF"}
    gzip_exts = {".gz", ".GZ", ".Gz"}

    if vcf_source.is_file():
        if vcf_source.suffix in gzip_exts:
            decompressed_path = vcf_source.with_suffix("")
        else:
            decompressed_path = vcf_source

        if decompressed_path.suffix not in vcf_exts:
            raise NotVCF(
                f"Input source should contain vcf extension: " f"{vcf_source}"
            )
        return [vcf_source]
    elif not vcf_source.exists():
        raise FileNotFoundError(vcf_source)
    elif not vcf_source.is_dir():
        raise ValueError(
            "Sources input must be a file of vcf.gz format, or a directory"
            f"containing vcf.gz files. {vcf_source} is invalid."
        )

    detected = []

    for extension in vcf_exts.union(
        {f"{x}{y}" for x in vcf_exts for y in gzip_exts}
    ):
        for vcf_file_path in vcf_source.glob(f"*{extension}"):
            detected.append(vcf_file_path)

    n_detected = len(detected)

    if n_detected > 0:
        print(f"Detected {n_detected} vcf files in {vcf_source.as_posix()}")
    else:
        raise NoVCFs(f"No VCF files detected in {vcf_source.as_posix()}")

    return detected


_RSCRIPTS_DIR = Directories.r_scripts()
_VCF_PARSER_R_PATH = _RSCRIPTS_DIR / "parse_vcf.R"


def annotate_source(
    source_path: Path,
    assembly: str,
    output_name: str | None = None,
    cache_directory: Path = Directories.app_annotated_inputs(),
):
    vcf_paths = find_vcf_files(source_path)
    target_filenames = [
        (v.with_suffix("").with_suffix(".vcf.anno")).name for v in vcf_paths
    ]
    target_paths = [cache_directory / tf for tf in target_filenames]

    for source, target in zip(vcf_paths, target_paths):
        pipe(
            [
                "Rscript",
                _VCF_PARSER_R_PATH.as_posix(),
                "-v",
                source.as_posix(),
                "-t",
                Directories.data().as_posix(),
                "-a",
                assembly,
                "-w",
                _RSCRIPTS_DIR.as_posix(),
                "-o",
                target.as_posix(),
            ]
        )

    output_path = cache_directory / f"{output_name}.vcf.anno"

    if output_name is not None:
        if len(vcf_paths) == 1:
            target_paths[0].rename(output_path)
        else:
            print(f"-- building merged file: {output_path}")

            with open(output_path, "w") as merged_file:
                for written_path in target_paths:
                    print(f"-> {written_path}")
                    with open(written_path, "r") as g:
                        lines = g.readlines()

                    if written_path != target_paths[-1]:
                        lines[-1] += "\n"

                    merged_file.writelines(lines)
