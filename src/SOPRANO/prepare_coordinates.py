import pathlib
import subprocess
import tempfile


def _filter_transcript_file(
    bed_file: pathlib.PosixPath,
    transcript_file: pathlib.PosixPath,
    cache_dir: pathlib.PosixPath,
) -> pathlib.PosixPath:
    """

    Implementation of methods in line 92-93

    cut -f1 $BED | sort -u |
        fgrep -w -f - $SUPA/ensemble_transcript_protein.length >
            $TMP/$NAME.protein_length_filt.txt

    :param bed_file: Bed file representation of target protein regions
    :param transcript_file: Transcript file to filter
    :param cache_dir: cache directory
    :return: PosixPath to filtered file
    """

    # Get transcript file name and build filtered filename
    transcript_file_name = transcript_file.name
    transcript_filt_name = transcript_file_name.replace(
        ".length", "_filt.length"
    )

    # Rebuild path for filtered file
    transcript_filt_file = cache_dir.joinpath(transcript_filt_name)

    # Create list of cmd line args and pass to subprocess
    cli_str_cmd = (
        f"cut -f1 {bed_file} | "
        f"sort -u | "
        f"fgrep -w -f - {transcript_file} > "
        f"{transcript_filt_file}"
    )
    cli_cmd = cli_str_cmd.split(sep=" ")
    subprocess.run(cli_cmd)

    return transcript_filt_file


def filter_transcript_files(
    bed_file: pathlib.PosixPath,
    ensembl_transcript_protein: pathlib.PosixPath,
    ensembl_transcript: pathlib.PosixPath,
):
    """

    Implementation of lines 92-93

    Get list of transcripts from annotated bed files, filtering out
    those transcripts not present in the database

    :param bed_file:
    :param ensembl_transcript_protein:
    :param ensembl_transcript:
    :return: Posix paths to filtered transcripts and temporary directory object
    """
    tmp = tempfile.TemporaryDirectory()
    tmp_dir = pathlib.PosixPath(tmp.name)

    ensemble_trans_prot_filt = _filter_transcript_file(
        bed_file, ensembl_transcript_protein, tmp_dir
    )
    ensemble_trans_filt = _filter_transcript_file(
        bed_file, ensembl_transcript, tmp_dir
    )

    return ensemble_trans_prot_filt, ensemble_trans_filt, tmp


def randomize_protein_positions(*args, **kwargs):
    """
    Implementation of line 96

    :param args:
    :param kwargs:
    :return:
    """
    pass


def randomize_target_regions(*args, **kwargs):
    """
    Implementatino of line 111

    :param args:
    :param kwargs:
    :return:
    """
    pass


def exclude_positively_selected_genes(*args, **kwargs):
    """
    Implmentation of line 128
    :param args:
    :param kwargs:
    :return:
    """
    pass


def get_protein_complement(*args, **kwargs):
    """
    Implmentatino of line 139
    :param args:
    :param kwargs:
    :return:
    """
    pass


def transform_protein_coordinates(*args, **kwargs):
    """
    Implementation of line 144
    :param args:
    :param kwargs:
    :return:
    """
    pass
