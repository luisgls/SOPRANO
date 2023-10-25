import argparse
import pathlib
from typing import List, Tuple

from SOPRANO.utils import sh_utils
from SOPRANO.utils.path_utils import Directories
from SOPRANO.utils.sh_utils import pipe


def _get_src_dst_link_pairs(vep_cache: pathlib.Path):
    """
    Identify vep sources that can be linked to the soprano cache

    vep (src) cache: /path/to/species/release/*
    soprano (dst) cache: /path/to/src/SOPRANO/data/species/release/*

    :param vep_cache: location of vep cache outside soprano
    :return: list of tuples pairing src and dst dirs
    """

    species_caches = [
        specie for specie in vep_cache.glob("*") if specie.is_dir()
    ]

    src_dirs = []
    for item in [
        [release for release in vs.glob("*") if release.is_dir()]
        for vs in species_caches
    ]:
        src_dirs += item

    _src_dirs = [
        pathlib.Path(vr.parts[-2]).joinpath(vr.parts[-1]) for vr in src_dirs
    ]

    dst_dirs = [Directories.data(sd) for sd in _src_dirs]

    src_dst_link_pairs = [
        (std_sys, soprano)
        for std_sys, soprano in zip(src_dirs, dst_dirs)
        if not soprano.exists()
    ]

    return src_dst_link_pairs


def _process(
    src_dir: pathlib.Path,
    dst_dir: pathlib.Path,
    pattern: str,
    compute_chroms: bool,
):
    fa_glob = list(src_dir.glob(f"*dna.{pattern}*.fa"))
    fai_glob = list(src_dir.glob(f"*dna.{pattern}*.fai"))

    src_fa_found = len(fa_glob) == 1
    src_fai_found = len(fai_glob) == 1

    if src_fa_found or src_fai_found:
        if not dst_dir.exists():
            print(f"Building directory: {dst_dir}")
            dst_dir.mkdir(parents=True)
    else:
        print(
            f"No appropriate fasta or fasta index files detected in {src_dir} "
            f"for {pattern}"
        )

    if src_fa_found:
        src_fa = fa_glob[0]
        dst_fa = dst_dir.joinpath(src_fa.name)

        if not dst_fa.exists():
            print(f"Sym linking {src_fa} -> {dst_fa}")
            dst_fa.symlink_to(src_fa)

        dst_chrom = dst_fa.with_suffix(".chrom")
        if not dst_chrom.exists() and compute_chroms:
            print(f"Computing chrom sizes: {dst_chrom}")
            pipe(["cut", "-f1,2", src_fa.as_posix()], output_path=dst_chrom)

    if src_fai_found:
        src_fai = fai_glob[0]
        dst_fai = dst_dir.joinpath(src_fai.name)
        print(f"Sym linking {src_fai} -> {dst_fai}")
        dst_fai.symlink_to(src_fai)


def _link_src_dst_pairs(
    src_dst_pairs: List[Tuple[pathlib.Path, pathlib.Path]],
    _skip_user_input=False,
):
    for src, dst in src_dst_pairs:
        str_response = "y" if _skip_user_input else ""

        input_options = ("y", "n")

        while str_response not in input_options:
            str_response = input(
                f"Do you wish to link {src} [y or n]: "
            ).lower()
            if str_response not in input_options:
                print(f"Invalid input: {str_response}\nEnter 'y' or 'n'\n")

        perform_link = str_response == "y"

        if perform_link:
            _process(src, dst, "toplevel", True)
            _process(src, dst, "primary_assembly", False)


def _link_vep_cache(vep_cache: pathlib.Path):
    if not vep_cache.exists():
        print(f"VEP cache not found: {vep_cache}")
    else:
        src_dst = _get_src_dst_link_pairs(vep_cache)
        _link_src_dst_pairs(src_dst)


def _link_vep_cache_parser():
    parser = argparse.ArgumentParser(description="VEP cache parser")

    parser.add_argument(
        "--cache",
        "-c",
        dest="src_cache",
        type=pathlib.Path,
        help="Provide the path to the ensembl vep cache. By default, will "
        "attempt to link sources from $HOME/.vep",
        default=Directories.std_sys_vep(),
    )

    src_cache: pathlib.Path = parser.parse_args().src_cache

    if not src_cache.exists():
        raise FileNotFoundError(src_cache)

    if not src_cache.is_dir():
        raise NotADirectoryError(src_cache)

    return src_cache


def annotate_homo_sapiens(
    input_file: pathlib.Path,
    assembly="GRCh38",
    output_file: pathlib.Path | None = None,
    overwrite=False,
    fasta_file=Directories.data("ensemble_transcriptID.translated.fasta"),
):
    # TODO: Complete and test ! This is place holding...

    if output_file is None:
        output_file = input_file.with_suffix(".anno")

    sh_utils.pipe(
        [
            "vep",
            "-i",
            input_file.as_posix(),
            "-o",
            output_file.as_posix(),
            "--cache",
            "--dir_cache",
            Directories.data().as_posix(),
            "--assembly",
            assembly,
            "--all_refseq",
            "--pick",
            "--symbol",
            "--no_stats",
            "--force_overwrite",  # TODO: Fix optional
            "--fasta",
            fasta_file.as_posix(),
        ]
    )


if __name__ == "__main__":
    # Check annotation procedure on vcf input
    test_input = Directories.examples("homo_sapiens_GRCh38.vcf")
    annotate_homo_sapiens(test_input)
