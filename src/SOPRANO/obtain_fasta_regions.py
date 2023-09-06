from SOPRANO import objects
from SOPRANO.pipeline_utils import MissingDataError, _PipelineComponent
from SOPRANO.sh_utils import subprocess_pipes


def _get_target_fasta_regions(
    paths: objects.AnalysisPaths, transcripts: objects.TranscriptPaths
):
    """
    Implements:

    bedtools getfasta -fi $TRANS -bed $TMP/$NAME.epitopes_cds.bed -fo
        $TMP/$NAME.epitopes_cds.fasta

    Extracts sequence from transcript ids fasta file for each of the
    intervals defined in the cds coords epitopes file.

    :param paths: AnalysisPaths instance, containing epitope cds path
                    definitions
    :param transcripts: TranscriptsPaths instance, containing ensembl
                    trancript id fasta file path definition
    """

    subprocess_pipes.pipe(
        [
            "bedtools",
            "getfasta",
            "-fi",
            transcripts.transcript_fasta.as_posix(),
            "-bed",
            paths.epitopes_cds.as_posix(),
            "-fo",
            paths.epitopes_cds_fasta,
        ]
    )


def _get_non_target_regions(
    paths: objects.AnalysisPaths, transcripts: objects.TranscriptPaths
):
    """
    Implements:

    bedtools getfasta -fi $TRANS -bed $TMP/$NAME.epitopes_cds.bed -fo
        $TMP/$NAME.epitopes_cds.fasta

    Extracts sequence from transcript ids fasta file for each of the
    intervals defined in the cds coords epitopes file.

    :param paths: AnalysisPaths instance, containing epitope cds path
                    definitions
    :param transcripts: TranscriptsPaths instance, containing ensembl
                    trancript id fasta file path definition
    """

    subprocess_pipes.pipe(
        [
            "bedtools",
            "getfasta",
            "-fi",
            transcripts.transcript_fasta.as_posix(),
            "-bed",
            paths.intra_epitopes_cds.as_posix(),
            "-fo",
            paths.epitopes_cds_fasta,
        ]
    )


class ObtainFastaRegions(_PipelineComponent):
    @staticmethod
    def check_ready(params: objects.Parameters):
        paths = (
            params.transcripts.transcript_fasta,
            params.epitopes_cds,
            params.intra_epitopes_cds,
        )

        for path in paths:
            if not path.exists():
                raise MissingDataError(path)

    @staticmethod
    def apply(params: objects.Parameters):
        _get_target_fasta_regions(params, params.transcripts)
        _get_non_target_regions(params, params.transcripts)


def get_transcript_regions_for_site_numbers(*args, **kwargs):
    """
    Implment lines 163-164
    :param args:
    :param kwargs:
    :return:
    """
    pass
