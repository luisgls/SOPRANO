import pathlib

from SOPRANO.objects import AnalysisPaths, TranscriptPaths
from SOPRANO.sh_utils import subprocess_pipes


def _get_target_fasta_regions(
    paths: AnalysisPaths, transcripts: TranscriptPaths
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
    paths: AnalysisPaths, transcripts: TranscriptPaths
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
            paths.intra_epitopes_cds_fasta,
        ]
    )


def _get_trans_regs(cds_fasta: pathlib.Path, output: pathlib.Path):
    """
    Implements:

    grep ">" $TMP/$NAME.epitopes_cds.fasta | sed 's/>//g' > $TMP/$NAME.listA

    List of transcript:trgions to estimate number of sites in epitopes

    :param cds_fasta: path to CDS coords fasta file
    :param output: path to output file
    """

    # subprocess_pipes.pipe(
    #     ["grep", '">"', cds_fasta.as_posix()],
    #     ["sed", '"s/>//g"'],
    #     output_path=output,
    # )

    # NOTE: We have passed the grep into the second stage of the pipe
    # since this caused issues in the subprocess (for some reason...)
    subprocess_pipes.pipe(
        ["sed", "s/>//g", cds_fasta.as_posix()],
        ["grep", ":"],
        output_path=output,
    )
