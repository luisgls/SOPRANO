import pathlib

from SOPRANO.objects import AnalysisPaths, Parameters
from SOPRANO.pipeline_utils import (
    MissingDataError,
    _PipelineComponent,
    is_empty,
)
from SOPRANO.sh_utils import subprocess_pipes


def _intersect_by_frequency(
    corrected_matrix: pathlib.Path, nans_output: pathlib.Path
):
    """
    Implements:

    awk '{NA+=$4}{NS+=$6}END{print NA"\t"NS}'
        $TMP/$NAME.final_corrected_matrix_A.txt > $TMP/$NAME.epitope_NaNs.txt

    :param corrected_matrix: Path to corrections for sites
        (epitopes of intra epitopes)
    :param nans_output: path to output for detected nans
    """

    subprocess_pipes.pipe(
        [
            "awk",
            r'{NA+=$4}{NS+=$6}END{print NA"\t"NS}',
            corrected_matrix.as_posix(),
        ],
        output_path=nans_output,
    )


class IntersectByFrequency(_PipelineComponent):
    @staticmethod
    def check_ready(params: Parameters):
        paths = (
            params.final_epitope_corrections,
            params.final_intra_epitope_corrections,
        )

        for path in paths:
            if not path.exists():
                raise MissingDataError(path)

    @staticmethod
    def apply(params: Parameters):
        IntersectByFrequency.check_ready(params)
        _intersect_by_frequency(
            params.final_epitope_corrections, params.epitope_nans
        )
        _intersect_by_frequency(
            params.final_intra_epitope_corrections, params.intra_epitope_nans
        )


_VARIANT_TYPES = (
    r"#|intergenic_variant|UTR|downstream|intron|miRNA|frameshift|"
    r"non_coding|splice_acceptor_variant|splice_donor_variant|"
    r"TF_binding_site_variant|upstream|incomplete|"
    r"regulatory_region_variant|retained|\?"
)

_FMT_COMMANDS = (
    ["awk", r'{if(length($3)>1||$10=="-"){}else{print}}'],
    ["cut", "-f4,5,7,10,89", "-"],
    ["sed", r"s/\//\t/g"],
    ["awk", r'{print $2"\t"$4"\t"$4"\t"$3}'],
    ["egrep", "-v", "-w", "-e", "coding_sequence_variant"],
    ["grep", "-v", "ENSEMBLTRANSCRIPT"],
)


def _get_silent_variant_counts(paths: AnalysisPaths):
    r"""
    Implements:

    egrep -v -e '#|intergenic_variant|UTR|downstream|intron|miRNA|frameshift|
        non_coding|splice_acceptor_variant|splice_donor_variant|
        TF_binding_site_variant|upstream|incomplete|regulatory_region_variant|
        retained|\?'
            $FILE | grep -w synonymous_variant |
                awk '{if(length($3)>1||$10=="-"){}else{print}}' |
                    cut -f4,5,7,10,89 -  | sed 's/\//\t/g' |
                        awk '{print $2"\t"$4"\t"$4"\t"$3}' |
                            egrep -v -w -e "coding_sequence_variant" |
                                grep -v "ENSEMBLTRANSCRIPT" >
                                    $TMP/$NAME.silent.bed

    :param paths:
    """

    subprocess_pipes.pipe(
        [
            "egrep",
            "-v",
            "-e",
            _VARIANT_TYPES,
            paths.input_path.as_posix(),
        ],
        ["grep", "-w", "synonymous_variant"],
        *_FMT_COMMANDS,
        output_path=paths.variants_silent,
    )


class GetSilentCounts(_PipelineComponent):
    @staticmethod
    def apply(params: Parameters):
        _get_silent_variant_counts(params)


def _get_nonsilent_variant_counts(paths: AnalysisPaths):
    r"""
    Implements:

    egrep -v -e '#|intergenic_variant|UTR|downstream|intron|miRNA|frameshift|
        non_coding|splice_acceptor_variant|splice_donor_variant|
        TF_binding_site_variant|upstream|incomplete|
        regulatory_region_variant|retained|\?' $FILE |
            grep -w -v synonymous_variant |
                awk '{if(length($3)>1||$10=="-"){}else{print}}'  |
                    cut -f4,5,7,10,89 -  |
                        sed 's/\//\t/g' |
                            awk '{print $2"\t"$4"\t"$4"\t"$3}' |
                                egrep -v -w -e "coding_sequence_variant" |
                                    grep -v "ENSEMBLTRANSCRIPT" >
                                        $TMP/$NAME.nonsilent.bed

    :param paths:
    """

    subprocess_pipes.pipe(
        [
            "egrep",
            "-v",
            "-e",
            _VARIANT_TYPES,
            paths.input_path.as_posix(),
        ],
        ["grep", "-w", "-v", "synonymous_variant"],
        *_FMT_COMMANDS,
        output_path=paths.variants_nonsilent,
    )


class GetNonSilentCounts(_PipelineComponent):
    @staticmethod
    def apply(params: Parameters):
        _get_nonsilent_variant_counts(params)


def _get_missense_variant_counts(paths: AnalysisPaths):
    r"""
    Implements:

    egrep -v -e '#|intergenic_variant|UTR|downstream|intron|miRNA|frameshift|
        non_coding|splice_acceptor_variant|splice_donor_variant|
        TF_binding_site_variant|upstream|incomplete|regulatory_region_variant|
        retained|\?' $FILE |
            grep -w -v synonymous_variant |
                grep -w missense_variant |
                    awk '{if(length($3)>1||$10=="-"){}else{print}}' |
                        cut -f4,5,7,10,89 -  |
                            sed 's/\//\t/g' |
                                awk '{print $2"\t"$4"\t"$4"\t"$3}' |
                                    egrep -v -w -e "coding_sequence_variant" |
                                        grep -v "ENSEMBLTRANSCRIPT" >
                                            $TMP/$NAME.missense.bed

    :param paths:
    :return:
    """
    subprocess_pipes.pipe(
        [
            "egrep",
            "-v",
            "-e",
            _VARIANT_TYPES,
            paths.input_path.as_posix(),
        ],
        ["grep", "-w", "-v", "synonymous_variant"],
        ["grep", "-w", "missense_variant"],
        *_FMT_COMMANDS,
        output_path=paths.variants_missense,
    )


class GetMissenseCounts(_PipelineComponent):
    @staticmethod
    def apply(params: Parameters):
        _get_missense_variant_counts(params)


def _get_intronic_variant_counts(paths: AnalysisPaths):
    """
    Implements:

    grep -v "^#" $NAME | grep -w "intron_variant" | grep -v "splice" |
        awk -F"\t|_" '{FS="\t|_"}{print $1"_"$7"\t"$2"\t"$2"\t"$3}' >
            $TMP/$NAME.intronic.bed

    NOTE: We added another step

    tr "_" "\t"

    to the implemented pipe before calling awk since the seperator was failing.

    :param paths:
    :return:
    """

    subprocess_pipes.pipe(
        ["grep", "-v", r"^#", paths.input_path.as_posix()],
        ["grep", "-w", "intron_variant"],
        ["grep", "-v", "splice"],
        ["tr", "_", "\t"],
        ["awk", '{print $1"_"$7"\t"$2"\t"$2"\t"$3}'],
        output_path=paths.variants_intronic,
    )


class GetIntronicCounts(_PipelineComponent):
    @staticmethod
    def apply(params: Parameters):
        _get_intronic_variant_counts(params)


def _count_mutations(variant_counts: pathlib.Path, output_path: pathlib.Path):
    if is_empty(variant_counts):
        counts = "0"
    else:
        counts = subprocess_pipes.pipe(
            ["wc", "-l", variant_counts.as_posix()],
            ["awk", r"{print $1}"],
            output_path=output_path,
        )

    return counts


def _count_intersected_mutations(
    variant_counts: pathlib.Path,
    on_off_regions: pathlib.Path,
    output_path: pathlib.Path,
):
    if is_empty(variant_counts):
        counts = "0"
    else:
        counts = subprocess_pipes.pipe(
            [
                "intersectBed",
                "-b",
                variant_counts.as_posix(),
                "-a",
                on_off_regions.as_posix(),
                "-wo",
            ],
            ["wc", "-l"],
            ["awk", r"{ print $1 }"],
            output_path=output_path,
        )

    return counts


class OnOffCounts(_PipelineComponent):
    @staticmethod
    def check_ready(params: Parameters):
        paths = (
            params.variants_silent,
            params.variants_nonsilent,
            params.variants_missense,
            params.variants_intronic,
            params.epitopes,
            params.intra_epitopes_prot,
        )

        for p in paths:
            if not p.exists():
                raise MissingDataError(p)

    @staticmethod
    def apply(params: Parameters):
        raw_silent = _count_mutations(
            params.variants_silent, params.raw_silent_count
        )
        raw_nonsilent = _count_mutations(
            params.variants_nonsilent, params.raw_nonsilent_count
        )
        raw_missense = _count_mutations(
            params.variants_missense, params.raw_missense_count
        )
        in_silent = _count_intersected_mutations(
            params.variants_silent, params.epitopes, params.in_silent_count
        )
        in_nonsilent = _count_intersected_mutations(
            params.variants_nonsilent,
            params.epitopes,
            params.in_nonsilent_count,
        )
        in_missense = _count_intersected_mutations(
            params.variants_missense, params.epitopes, params.in_missense_count
        )
        out_silent = _count_intersected_mutations(
            params.variants_silent,
            params.intra_epitopes_prot,
            params.out_silent_count,
        )
        out_nonsilent = _count_intersected_mutations(
            params.variants_nonsilent,
            params.intra_epitopes_prot,
            params.out_nonsilent_count,
        )
        out_missense = _count_intersected_mutations(
            params.variants_missense,
            params.intra_epitopes_prot,
            params.out_missense_count,
        )

        def _print(region, silent, nonsilent, missense):
            print(f"{region}:")
            print("{0:.<30}".format("Silent") + silent)
            print("{0:.<30}".format("Non-silent") + nonsilent)
            print("{0:.<30}".format("Missense") + missense)

        _print("Global region", raw_silent, raw_nonsilent, raw_missense)
        _print("(ON) Target region", in_silent, in_nonsilent, in_missense)
        _print("(OFF) Target region", out_silent, out_nonsilent, out_missense)
