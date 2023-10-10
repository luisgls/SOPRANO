import pathlib

from SOPRANO.utils.path_utils import Directories
from SOPRANO.utils.sh_utils import pipe


def join_hla_alleles(*hla_alleles):
    updated_hla_alleles = []

    for allele in hla_alleles:
        if not allele.startswith("HLA-"):
            allele = f"HLA-{allele}".upper()
        updated_hla_alleles.append(allele)

    return "|".join(updated_hla_alleles) + "|"


def immunopeptidome_from_hla(
    *hla_alleles: str, output_name: str, cache_loc: pathlib.Path | None = None
):
    joined_alleles = join_hla_alleles(*hla_alleles)

    if not output_name.endswith(".bed"):
        output_name = f"{output_name}.bed"

    if cache_loc is None:
        cache_loc = Directories.app_immunopeptidomes()

    output_path = cache_loc / output_name

    hla_binders_path = Directories.data("allhlaBinders_exprmean1.IEDBpeps.bed")

    pipe(
        ["egrep", "-w", "-e", joined_alleles, hla_binders_path.as_posix()],
        ["sortBed", "-i", "stdin"],
        ["mergeBed", "-i", "stdin"],
        output_path=output_path,
    )
