import pathlib
from typing import List

from SOPRANO.utils.path_utils import Directories
from SOPRANO.utils.sh_utils import pipe


def join_hla_alleles(*hla_alleles):
    updated_hla_alleles = []

    for allele in hla_alleles:
        if not allele.startswith("HLA-"):
            allele = f"HLA-{allele}".upper()
        updated_hla_alleles.append(allele)

    return r"\|".join(updated_hla_alleles)


def join_transcript_ids(transcript_ids: List[str]):
    return r"\|".join(transcript_ids)


def prior_filter_restrictions(
    transcript_ids: List[str], output_path: pathlib.Path
):
    tmp_path = output_path.with_suffix(".tmp")
    hla_binders_path = Directories.immunopeptidome_aux_files(
        "allhlaBinders_exprmean1.IEDBpeps.bed"
    )
    print(f"Retaining transcripts {join_transcript_ids(transcript_ids)}")
    pipe(
        [
            "grep",
            "-w",
            "-e",
            join_transcript_ids(transcript_ids),
            hla_binders_path.as_posix(),
        ],
        output_path=tmp_path,
        overwrite=True,
    )
    return tmp_path


def prior_filter_exclusions(
    transcript_ids: List[str], output_path: pathlib.Path | None = None
):
    assert output_path is not None

    tmp_path = output_path.with_suffix(".tmp")
    # TODO: Should this be a runtime option?
    hla_binders_path = Directories.immunopeptidome_aux_files(
        "allhlaBinders_exprmean1.IEDBpeps.bed"
    )
    print(f"Excluding transcripts {join_transcript_ids(transcript_ids)}")
    pipe(
        [
            "grep",
            "-v",
            "-w",
            "-e",
            join_transcript_ids(transcript_ids),
            hla_binders_path.as_posix(),
        ],
        output_path=tmp_path,
        overwrite=True,
    )
    return tmp_path


def immunopeptidome_from_hla(
    *hla_alleles: str,
    output_name: str,
    cache_loc: pathlib.Path | None = None,
    restricted_transcript_ids: List[str] = [],
    excluded_transcript_ids: List[str] = [],
):
    if not output_name.endswith(".bed"):
        output_name = f"{output_name}.bed"

    if cache_loc is None:
        cache_loc = Directories.app_immunopeptidomes()

    output_path = cache_loc / output_name

    n_restricted = len(restricted_transcript_ids)
    n_excluded = len(excluded_transcript_ids)

    if (n_restricted > 0) and (n_excluded > 0):
        raise ValueError(
            "Cannot restrict and exclude transcripts simultaneously"
        )
    elif n_restricted > 0:
        use_input = prior_filter_restrictions(
            restricted_transcript_ids, output_path=output_path
        )
    elif n_excluded > 0:
        use_input = prior_filter_exclusions(
            excluded_transcript_ids, output_path=output_path
        )
    else:
        # TODO: Should this be an option, e.g. use netMHCpan version, etc
        use_input = Directories.immunopeptidome_aux_files(
            "allhlaBinders_exprmean1.IEDBpeps.bed"
        )

    joined_alleles = join_hla_alleles(*hla_alleles)

    print(f"Filtering by alleles: {joined_alleles}")

    pipe(
        ["grep", "-w", "-e", joined_alleles, use_input.as_posix()],
        ["sortBed", "-i", "stdin"],
        ["mergeBed", "-i", "stdin"],
        output_path=output_path,
    )

    tmp_path = output_path.with_suffix(".tmp")
    tmp_path.unlink(missing_ok=True)
    print(f"All done: {output_path.as_posix()}")
