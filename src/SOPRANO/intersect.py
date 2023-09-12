import pathlib

from SOPRANO.objects import AnalysisPaths, Parameters
from SOPRANO.pipeline_utils import MissingDataError, _PipelineComponent
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


def _get_silent_variant_counts(paths: AnalysisPaths):
    """
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
            r"#|intergenic_variant|UTR|downstream|intron|miRNA|frameshift|"
            r"non_coding|splice_acceptor_variant|splice_donor_variant|"
            r"TF_binding_site_variant|upstream|incomplete|"
            r"regulatory_region_variant|retained|\?",
            paths.input_path.as_posix(),
        ],
        ["grep", "-w", "synonymous_variant"],
        ["awk", r'{if(length($3)>1||$10=="-"){}else{print}}'],
        ["cut", "-f4,5,7,10,89", "-"],
        ["sed", r"s/\//\t/g"],
        ["awk", r'{print $2"\t"$4"\t"$4"\t"$3}'],
        ["egrep", "-v", "-w", "-e", "coding_sequence_variant"],
        ["grep", "-v", "ENSEMBLTRANSCRIPT"],
        output_path=paths.variants_silent,
    )


class GetSilentCounts(_PipelineComponent):
    @staticmethod
    def apply(params: Parameters):
        _get_silent_variant_counts(params)
