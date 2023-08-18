import pathlib
import tempfile


def filter_bed_file_transcript(
    bed_file: pathlib.PosixPath,
    transcript_file: pathlib.PosixPath,
    cache_dir: pathlib.PosixPath,
):
    """
    cut -f1 $BED | sort -u |
        fgrep -w -f - $SUPA/ensemble_transcript_protein.length >
            $TMP/$NAME.protein_length_filt.txt

    :param bed_file: Bed file representation of target protein regions
    :param transcript_file: Transcript file to filter
    :param cache_dir: cache directory
    :return: PosixPath to filtered file
    """
    pass


def prepare_coordinate_files(
    bed_file: pathlib.PosixPath,
    ensembl_transcript_protein: pathlib.PosixPath,
    ensembl_transcript: pathlib.PosixPath,
):
    """
    Get list of transcripts from annotated bed files, filtering out
    those transcripts not present in the database

    :param bed_file:
    :param ensembl_transcript_protein:
    :param ensembl_transcript:
    :return: Posix paths to filtered transcripts and temporary directory object
    """
    tmp = tempfile.TemporaryDirectory()
    tmp_dir = pathlib.PosixPath(tmp.name)

    ensemble_trans_prot_filt = filter_bed_file_transcript(
        bed_file, ensembl_transcript_protein, tmp_dir
    )
    ensemble_trans_filt = filter_bed_file_transcript(
        bed_file, ensembl_transcript, tmp_dir
    )

    return ensemble_trans_prot_filt, ensemble_trans_filt, tmp
