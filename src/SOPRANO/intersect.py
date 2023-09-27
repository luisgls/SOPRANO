import pathlib

from SOPRANO.objects import AnalysisPaths, Parameters
from SOPRANO.pipeline_utils import (
    MissingDataError,
    SOPRANOError,
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
            ["awk", r"{print $1+1}"],
            output_path=output_path,
        )
        # In the awk print we have $1+1 since our internal piping
        # convention is such that when pipes are written to a file,
        # they do not have a trailing \n.
        # This means that wc -l will in fact count n_lines - 1
        # since it looks for \n

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


def get_counts(in_out_counts: pathlib.Path):
    return int(subprocess_pipes.pipe(["cat", in_out_counts.as_posix()]))


def _update_epitopes_data_file(
    variant_counts: pathlib.Path,
    in_out_counts: pathlib.Path,
    paths: AnalysisPaths,
    _use_epitope: bool,
    _label: str,
):
    """
    Implements:

    intersectBed -b $TMP/$NAME.nonsilent.bed -a $TMP/$NAME.epitopes.bed -wo |
        awk '{OFS="\t"}{print $1,"1","2",$0}' | cut -f1-11 -|
            sortBed -i stdin | mergeBed -i stdin -c 11 -o count | cut -f1,4 |
                awk '{print $0"\textra_missense_variant"}' >>
                    $TMP/$NAME.data_epitopes

    and

    intersectBed -b $TMP/$NAME.silent.bed -a
    $TMP/$NAME.intra_epitopes_prot.bed -wo |
        awk '{OFS="\t"}{print $1,"1","2",$0}' |
            cut -f1-11 -| sortBed -i stdin | mergeBed -i stdin -c 10 -o count |
                cut -f1,4 |
                    awk '{print $0"\tintra_synonymous_variant"}' >>
                        $TMP/$NAME.data_epitopes

    :param variant_counts:
    :param in_out_counts:
    :param paths:
    :return:
    """

    if get_counts(in_out_counts) > 0:
        if _use_epitope:
            epi_path = paths.epitopes
            col_idx = "11"
            awk_str = r'{print $0"\textra_' + _label + '_variant"}'
        else:
            epi_path = paths.intra_epitopes_prot
            col_idx = "10"
            awk_str = r'{print $0"\tintra_' + _label + '_variant"}'

        subprocess_pipes.pipe(
            [
                "intersectBed",
                "-b",
                variant_counts.as_posix(),
                "-a",
                epi_path.as_posix(),
                "-wo",
            ],
            ["awk", '{OFS="\t"}{print $1,"1","2",$0}'],
            ["cut", "-f1-11", "-"],
            ["sortBed", "-i", "stdin"],
            [
                "mergeBed",
                "-i",
                "stdin",
                "-c",
                col_idx,
                "-o",
                "count",
            ],
            ["cut", "-f1,4"],
            ["awk", awk_str],
            output_path=paths.data_epitopes,
            mode="a",
            overwrite=True,
        )


class BuildEpitopesDataFile(_PipelineComponent):
    @staticmethod
    def check_ready(params: Parameters):
        paths = (
            params.variants_silent,
            params.variants_nonsilent,
            params.in_silent_count,
            params.in_nonsilent_count,
            params.out_silent_count,
            params.out_nonsilent_count,
        )
        for path in paths:
            if not path.exists():
                raise MissingDataError(path)

    @staticmethod
    def apply(params: Parameters):
        BuildEpitopesDataFile.check_ready(params)

        variants = (params.variants_silent, params.variants_nonsilent)
        in_counts = (params.in_silent_count, params.in_nonsilent_count)
        out_counts = (params.out_silent_count, params.out_nonsilent_count)
        labs = ("synonymous", "missense")

        for variant_count, in_out_count, use_epi, lab in zip(
            variants * 2,
            in_counts + out_counts,
            (True, True, False, False),
            labs * 2,
        ):
            _update_epitopes_data_file(
                variant_count,
                in_out_count,
                params,
                _use_epitope=use_epi,
                _label=lab,
            )


def _check_target_mutations(paths: AnalysisPaths):
    in_silent_count = get_counts(paths.in_silent_count)
    in_nonsilent_count = get_counts(paths.in_nonsilent_count)
    in_missense_count = get_counts(paths.in_missense_count)

    if in_silent_count + in_nonsilent_count + in_missense_count == 0:
        raise SOPRANOError(
            f"No mutations found in target region for input file "
            f"{paths.input_path}"
        )


class CheckTargetMutations(_PipelineComponent):
    @staticmethod
    def check_ready(params: Parameters):
        paths = (
            params.in_silent_count,
            params.in_nonsilent_count,
            params.in_missense_count,
        )

        for path in paths:
            if not path.exists():
                raise MissingDataError(path)

    @staticmethod
    def apply(params: Parameters):
        CheckTargetMutations.check_ready(params)
        _check_target_mutations(params)
