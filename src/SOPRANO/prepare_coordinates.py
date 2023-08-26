import pathlib

from SOPRANO.sh_utils import subprocess_pipes


def _filter_transcript_file(
    bed_file: pathlib.Path,
    transcript_file: pathlib.Path,
    cache_dir: pathlib.Path,
) -> pathlib.Path:
    """

    Implementation of methods in line 92-93

    cut -f1 $BED | sort -u |
        fgrep -w -f - $SUPA/ensemble_transcript_protein.length >
            $TMP/$NAME.protein_length_filt.txt

    Bedfile format looks like:
    ENST00000000233 9       18
    ENST00000000233 37      46
    ENST00000000233 98      111
    ENST00000000233 115     124
    ENST00000000233 164     177
    ...

    cut -f1 $BED extracts first col, e.g.
    ENST00000000233
    ENST00000000233
    ENST00000000233
    ENST00000000233
    ENST00000000233

    sort -u extracts the unique lines

    fgrep -w -f <file>, -w matches only whole words, -f take patterns from file

    Hence, the operation strips the first column ENST<number>, looks for all
    the unique identifiers, then regex matches the identifiers against those
    in the transcript files.

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

    # Perform filtering
    subprocess_pipes.pipe(
        ["cut", "f1", bed_file.as_posix()],
        ["sort", "-u"],
        ["fgrep", "-w", "-f", "-", transcript_file.as_posix()],
        output_path=transcript_filt_file,
    )

    return transcript_filt_file


def filter_transcript_files(
    bed_file: pathlib.Path,
    ensembl_transcript_protein: pathlib.Path,
    ensembl_transcript: pathlib.Path,
    cache_dir: pathlib.Path,
):
    """

    Implementation of lines 92-93

    Get list of transcripts from annotated bed files, filtering out
    those transcripts not present in the database

    :param bed_file: input bedfile
    :param ensembl_transcript_protein: ensembl transcript for protein
    :param ensembl_transcript: ensembl transcript for (?)
    :param cache_dir: directory to cache results into
    :return: Posix paths to filtered transcripts
    """

    ensemble_trans_prot_filt = _filter_transcript_file(
        bed_file, ensembl_transcript_protein, cache_dir
    )
    ensemble_trans_filt = _filter_transcript_file(
        bed_file, ensembl_transcript, cache_dir
    )

    return ensemble_trans_prot_filt, ensemble_trans_filt


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
