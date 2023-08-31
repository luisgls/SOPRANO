import pathlib
import subprocess

from SOPRANO.objects import AnalysisPaths, AuxiliaryPaths, TranscriptPaths
from SOPRANO.sh_utils import subprocess_pipes


class DependentDataError(Exception):
    pass


def _filter_transcript_file(
    bed_file: pathlib.Path,
    transcript_file: pathlib.Path,
    transcript_filt: pathlib.Path,
) -> None:
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
    in the transcript files. (we are finding the matches IN the transcript)

    :param bed_file: Bed file representation of target protein regions
    :param transcript_file: Transcript file to filter
    :param cache_dir: cache directory
    :return: PosixPath to filtered file
    """

    # Perform filtering
    subprocess_pipes.pipe(
        ["cut", "-f1", bed_file.as_posix()],
        ["sort", "-u"],
        ["fgrep", "-w", "-f", "-", transcript_file.as_posix()],
        output_path=transcript_filt,
    )


def filter_transcript_files(
    path: AnalysisPaths, transcripts: TranscriptPaths
) -> None:
    """
    Implementation of lines 92-93

    Get list of transcripts from annotated bed files, filtering out
    those transcripts not present in the database

    :param path: AnalysisPaths instance (contains e.g. bedfile path)
    :param transcripts: TranscriptsPath instance
    :return:
    """
    bed_file = path.bed_path
    transcript_protein_path = transcripts.protein_transcript_length
    transcript_path = transcripts.transcript_length

    transcript_protein_filt = path.filtered_protein_transcript
    transcript_filt = path.filtered_transcript

    _filter_transcript_file(
        bed_file,
        transcript_protein_path,
        transcript_protein_filt,
    )

    _filter_transcript_file(
        bed_file,
        transcript_path,
        transcript_filt,
    )


def _define_excluded_regions_for_randomization(
    paths: AnalysisPaths,
):
    """
    We want to execute the commands
    cut -f1,2,3 $BED > $TMP/$NAME.exclusion.ori
    cut -f1 $BED |
        awk '{OFS="\t"}{print $1,0,2}' |
            sortBed -i stdin  >> $TMP/$NAME.exclusion.ori

    the first command simply extracts the first 3 cols of the bed file
    and stores them in a tmp file *.exclusion.ori

    the second command pipes the following:
    1 - cut -f1 $BED
        extract the first col from bed file; then
    2 - awk '{OFS="\t"}{print $1,0,2}'
        uses awk to process the output from 1;
        {OFS='\t'} tab delimits the output
        {print $1,0,2} prints the first column ($1) then
        delimits the numbers 0 and 2

        Summary: this pipe tab delimits the first col, followed
        by the numbers 0 and 2
    3 - sortBed -i stdin
        sortBed sorts input files by features (e.g. chrom size).
        by passing -i stdin, this is telling bed to use the
        standard input stream. Therefore, it can be used in pipes.

        By default, sorts a BED file by chromosome and then by start position
        in ascending order.

    Net result of this function:

    Take a bed file and extract the first 3 cols.
    Append to the bottom of this file the sorted list of chromosomes,
    followed by the numbers 0 and 2.

    :param name:
    :param bed_path:
    :param tmpdir:
    :return:
    """

    cut_bed_proc = subprocess.run(
        ["cut", "-f1,2,3", paths.bed_path.as_posix()], capture_output=True
    )

    subprocess_pipes.process_output_to_file(
        cut_bed_proc, path=paths.exclusions
    )

    subprocess_pipes.pipe(
        ["cut", "-f1", paths.bed_path.as_posix()],
        ["awk", '{OFS="\t"}{print $1,0,2}'],
        ["sortBed", "-i", "stdin"],
        output_path=paths.exclusions,
        mode="a",
        overwrite=True,
    )


def _sort_excluded_regions_for_randomization(
    paths: AnalysisPaths, seed: int | None = None
) -> None:
    """
    Implement
    sortBed -i $TMP/$NAME.exclusion.ori > $TMP/$NAME.exclusion.bed
    bedtools shuffle -i $BED -g $TMP/$NAME.protein_length_filt.txt
        -excl $TMP/$NAME.exclusion.bed -chrom > $TMP/$NAME.epitopes.ori2


    :param paths: AnalaysisPaths instance
    """

    # Sort list of exclusions obtained from
    # _define_excluded_regions_for_randomization
    subprocess_pipes.pipe(
        ["sortBed", "-i", paths.exclusions.as_posix()],
        output_path=paths.exclusions_sorted,
    )

    # randomly permute the genomic locations of a feature file among a genome
    # defined in a genome file. Here the feature is the input bed file,
    # and the genome is that of the filtered protein file.
    # The -excl flag indicates to bedtools that we want to exclude regions
    # from the shuffling procedure.
    # The -chrom flag tells bedtools to keep features in -i on the same
    # chromosome. Solely permute their location on the chromosome.
    pipe_args = [
        "bedtools",
        "shuffle",
        "-i",
        paths.bed_path.as_posix(),
        "-g",
        paths.filtered_protein_transcript,
        "-excl",
        paths.exclusions_sorted,
        "-chrom",
    ]

    if seed is not None:
        pipe_args += ["-seed", str(seed)]

    subprocess_pipes.pipe(
        pipe_args,
        output_path=paths.exclusions_shuffled,
    )


def _randomize_with_target_file(
    paths: AnalysisPaths, transcripts: TranscriptPaths, seed: int | None = None
):
    """
    Implementation of
    cut -f1 $TARGET | sort -u | fgrep -w -f -
        $SUPA/ensemble_transcript_protein.length >>
            $TMP/$NAME.protein_length_filt.txt

    cut -f1 $TARGET | sort -u | fgrep -w -f -
        $SUPA/ensemble_transcript.length >>
            $TMP/$NAME.transcript_length_filt.txt

    These two steps take the target region file, extract the unique chroms,
    find those which exist in the transcript file(s), then append these data
    to the end of the filtered transcript file(s)

    bedtools shuffle -i $BED -g $TMP/$NAME.protein_length_filt.txt -incl
        $TARGET -noOverlapping > $TMP/$NAME.epitopes.ori2

    This ensures that that all shuffling takes place across the TARGET file,
    and that no intervals occupy a single common base pair

    :param paths:
    :param seed: rng seed
    :return:
    """

    if paths.target_regions_path is None:
        raise ValueError("Method only valid if target regions is not none!")

    subprocess_pipes.pipe(
        ["cut", "-f1", paths.target_regions_path.as_posix()],
        ["sort", "-u"],
        [
            "fgrep",
            "-w",
            "-f",
            "-",
            transcripts.protein_transcript_length.as_posix(),
        ],
        output_path=paths.filtered_protein_transcript,
        mode="a",
        overwrite=True,
    )

    subprocess_pipes.pipe(
        ["cut", "-f1", paths.target_regions_path.as_posix()],
        ["sort", "-u"],
        ["fgrep", "-w", "-f", "-", transcripts.transcript_length.as_posix()],
        output_path=paths.filtered_transcript,
        mode="a",
        overwrite=True,
    )

    pipe_args = [
        "bedtools",
        "shuffle",
        "-i",
        paths.bed_path.as_posix(),
        "-g",
        paths.filtered_protein_transcript.as_posix(),
        "-incl",
        paths.target_regions_path.as_posix(),
        "-noOverlapping",
    ]

    if seed is not None:
        pipe_args += ["-seed", str(seed)]

    subprocess_pipes.pipe(
        pipe_args,
        output_path=paths.exclusions_shuffled,
    )

    # TODO: Finish unit testing... This is quite tricky


def _non_randomized(paths: AnalysisPaths):
    """

    TODO: Shouldn't this be instead

    ----------

    sort -k 1,1 -k2,2n -u $BED > ...

    Reasoning: Consider the following case

    ENST00000001008 189     198
    ENST00000001008 27      36

    If sort -u then this is already sorted!
    (presuming heirarchy of sorting should be chrom, then start, then stop)

    Applying sort -k 1,1 -k2,2n -u $BED > ...

    ENST00000001008 27      36
    ENST00000001008 189     198

    ----------

    Implement
    sort -u $BED > $TMP/$NAME.epitopes.ori2
    :param paths:
    :return:
    """
    subprocess_pipes.pipe(
        ["sort", "-u", paths.bed_path], output_path=paths.exclusions_shuffled
    )


def randomize_protein_positions(*args, **kwargs):
    """
    Implementation of line 96

    :param args:
    :param kwargs:
    :return:
    """
    pass


def _exclude_positively_selected_genes_disabled(paths: AnalysisPaths):
    """
    Implement
    cp $TMP/$NAME.epitopes.ori2 $TMP/$NAME.epitopes.bed
    :param paths:
    :return:
    """
    subprocess_pipes.pipe(["cp", paths.exclusions_shuffled, paths.epitopes])


def _exclude_positively_selected_genes(paths: AnalysisPaths):
    """
    Implement
    fgrep -w -v -f $SUPA/genes2exclude.txt $TMP/$NAME.epitopes.ori2 >
        $TMP/$NAME.epitopes.bed
    :param paths:
    :return:
    """

    subprocess_pipes.pipe(
        [
            "fgrep",
            "-w",
            "-v",
            "-f",
            AuxiliaryPaths.genes_to_exclude.as_posix(),
            paths.exclusions_shuffled.as_posix(),
        ],
        output_path=paths.epitopes,
    )

    # TODO: Unit test


def get_protein_complement(paths: AnalysisPaths):
    """
    Implement
    sortBed -i $TMP/$NAME.epitopes.bed |
        complementBed -i stdin -g $TMP/$NAME.protein_length_filt.txt >
            $TMP/$NAME.intra_epitopes_prot.tmp

    cut -f1 $TMP/$NAME.epitopes.bed |
        sort -u |
            fgrep -w -f - $TMP/$NAME.intra_epitopes_prot.tmp >
                $TMP/$NAME.intra_epitopes_prot.bed

    :param paths
    :return:
    """

    temporary_file = paths.intra_epitopes_prot.with_suffix(".tmp")

    subprocess_pipes.pipe(
        ["sortBed", "-i", paths.epitopes],
        [
            "complementBed",
            "-i",
            "stdin",
            "-g",
            paths.filtered_protein_transcript,
        ],
        output_path=temporary_file,
    )

    subprocess_pipes.pipe(
        ["cut", "-f1", paths.epitopes],
        ["sort", "-u"],
        ["fgrep", "-w", "-f", "-", temporary_file],
        output_path=paths.intra_epitopes_prot,
    )

    # TODO: unit test


def _prep_ssb192(paths: AnalysisPaths):
    """
    Implement
    awk '{OFS="\t"}{if( (($2*3)-6) >= 0 )
        {print $1,($2*3)-6,$3*3+3,$0}
            else{print $1,($2*3)-3,$3*3+3,$0}}' $TMP/$NAME.epitopes.bed >
                $TMP/$NAME.epitopes_cds.bed
    :param paths:
    :return:
    """

    subprocess_pipes.pipe(
        [
            "awk",
            '{OFS="\t"}{if( (($2*3)-6) >= 0 ) '
            "{print $1,($2*3)-6,$3*3+3,$0} "
            "else{print $1,($2*3)-3,$3*3+3,$0}}",
            paths.epitopes.as_posix(),
        ],
        output_path=paths.epitopes_cds,
    )

    # TODO: Unit test


def _prep_not_ssb192(paths: AnalysisPaths):
    """
    Implement
    awk '{OFS="\t"}{print $1,($2*3)-3,$3*3,$0}' $TMP/$NAME.epitopes.bed >
        $TMP/$NAME.epitopes_cds.bed
    :param paths:
    :return:
    """

    subprocess_pipes.pipe(
        ["awk", '{OFS="\t"}{print $1,($2*3)-3,$3*3,$0}', paths.epitopes],
        output_path=paths.epitopes_cds,
    )

    # TODO: Unit test


def transform_protein_coordinates(paths: AnalysisPaths):
    """
    Implement

    sortBed -i $TMP/$NAME.epitopes_cds.bed |
        complementBed -i stdin -g $TMP/$NAME.transcript_length_filt.txt >
            $TMP/$NAME.intra_epitopes.tmp

    cut -f1 $TMP/$NAME.epitopes.bed |
        sort -u |
            fgrep -w -f - $TMP/$NAME.intra_epitopes.tmp >
                $TMP/$NAME.intra_epitopes_cds.bed

    :param paths:
    :return:
    """

    temporary_file = paths.intra_epitopes.with_suffix(".tmp")

    subprocess_pipes.pipe(
        ["sortBed", "-i", paths.epitopes_cds],
        ["complementBed", "-i", "stdin", "-g", paths.filtered_transcript],
        output_path=temporary_file,
    )

    subprocess_pipes.pipe(
        ["cut", "-f1", paths.epitopes],
        ["sort", "-u"],
        ["fgrep", "-w", "-f", "-", temporary_file],
        output_path=paths.intra_epitopes,
    )

    # TODO: Unit test
