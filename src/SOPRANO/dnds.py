from SOPRANO.objects import (
    AnalysisPaths,
    AuxiliaryFiles,
    AuxiliaryPaths,
    Parameters,
)
from SOPRANO.pipeline_utils import MissingDataError, _PipelineComponent
from SOPRANO.sh_utils import subprocess_pipes


def _intersect_introns(paths: AnalysisPaths, aux_files: AuxiliaryPaths):
    """
    Implements:

    intersectBed -a $BASEDIR/data/transcript_intron_length.bed -b
    $TMP/$NAME.intronic.bed -wo |
        mergeBed -i stdin -c 4,5,6,10,11 -o mode,mode,mode,collapse,count |
            awk '{print $4"\t"$8/($6+1)"\t"$8"\t"$6}' >
                $TMP/$NAME.intronic.rate

    :param paths:
    :return:
    """

    subprocess_pipes.pipe(
        [
            "intersectBed",
            "-a",
            aux_files.intron_length.as_posix(),
            "-b",
            paths.variants_intronic.as_posix(),
            "-wo",
        ],
        [
            "mergeBed",
            "-i",
            "stdin",
            "-c",
            "4,5,6,10,11",
            "-o",
            "mode,mode,mode,collapse,count",
        ],
        ["awk", r'{print $4"\t"$8/($6+1)"\t"$8"\t"$6}'],
        output_path=paths.intron_rate,
    )


class ComputeIntronRate(_PipelineComponent):
    @staticmethod
    def check_ready(params: Parameters):
        if not params.variants_intronic.exists():
            raise MissingDataError(params.variants_intronic)

    @staticmethod
    def apply(params: Parameters):
        _intersect_introns(params, AuxiliaryFiles)
