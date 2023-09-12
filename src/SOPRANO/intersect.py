import pathlib

from SOPRANO.objects import Parameters
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
