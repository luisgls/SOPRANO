import pathlib
import subprocess

from SOPRANO.objects import (
    AnalysisPaths,
    AuxiliaryFiles,
    AuxiliaryPaths,
    Parameters,
    TranscriptPaths,
)
from SOPRANO.pipeline_utils import (
    MissingDataError,
    _PipelineComponent,
    _Randomize,
)
from SOPRANO.sh_utils import subprocess_pipes


def _filter_transcript_file(
    bed_file: pathlib.Path,
    transcript_file: pathlib.Path,
    filtered_transcript_file: pathlib.Path,
) -> None:
    """
    Implements:

    cut -f1 $BED | sort -u |
        fgrep -w -f - $SUPA/ensemble_transcript_protein.length >
            $TMP/$NAME.protein_length_filt.txt

    The operation strips the first column ENST<number>, looks for all
    the unique identifiers, then matches the identifiers against those
    in the transcript files.

    :param bed_file: Bed file representation of target protein regions
    :param transcript_file: Path to input transcript file
    :param filtered_transcript_file: Path to output filtered transcript
    """

    # Perform filtering
    subprocess_pipes.pipe(
        ["cut", "-f1", bed_file.as_posix()],
        ["sort", "-u"],
        ["fgrep", "-w", "-f", "-", transcript_file.as_posix()],
        output_path=filtered_transcript_file,
    )


def filter_transcript_files(
    paths: AnalysisPaths, transcripts: TranscriptPaths
) -> None:
    """
    Get list of transcripts from annotated bed files, filtering out
    those transcripts not present in the database

    :param paths: AnalysisPaths instance (contains e.g. bed_path attribute)
    :param transcripts: TranscriptsPath instance
    """
    bed_file = paths.bed_path
    transcript_protein_path = transcripts.protein_transcript_length
    transcript_path = transcripts.transcript_length

    transcript_protein_filt = paths.filtered_protein_transcript
    transcript_filt = paths.filtered_transcript

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


def _define_excluded_regions_for_randomization(paths: AnalysisPaths):
    """
    Implements:

    cut -f1,2,3 $BED > $TMP/$NAME.exclusion.ori
    cut -f1 $BED |
        awk '{OFS="\t"}{print $1,0,2}' |
            sortBed -i stdin  >> $TMP/$NAME.exclusion.ori

    Take a bed file and extract the first 3 cols.
    Append to the bottom of this file the sorted list of chromosomes,
    followed by the numbers 0 and 2.

    :param paths: Analysis paths instance containing bed file definition
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
    Implements:

    sortBed -i $TMP/$NAME.exclusion.ori > $TMP/$NAME.exclusion.bed
    bedtools shuffle -i $BED -g $TMP/$NAME.protein_length_filt.txt
        -excl $TMP/$NAME.exclusion.bed -chrom > $TMP/$NAME.epitopes.ori2

    Randomly permute the genomic locations within the input bed file among
    the filtered protein length file. We apply pre-defined exclusions from
    the end-points of shuffling. By invoking "-chrom", shuffled features
    are restricted to a location on the same chromosome.

    :param paths: AnalaysisPaths instance
    :param seed: Random number generator seed
    """

    subprocess_pipes.pipe(
        ["sortBed", "-i", paths.exclusions.as_posix()],
        output_path=paths.exclusions_sorted,
    )

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
    Implements:

    cut -f1 $TARGET | sort -u | fgrep -w -f -
        $SUPA/ensemble_transcript_protein.length >>
            $TMP/$NAME.protein_length_filt.txt

    cut -f1 $TARGET | sort -u | fgrep -w -f -
        $SUPA/ensemble_transcript.length >>
            $TMP/$NAME.transcript_length_filt.txt

    Takes the target region file, extracts the unique chroms, finds those which
    exist in the transcript file(s), then append these data to the end of the
    filtered transcript file(s).

    All shuffling takes place across the TARGET file, such that no intervals
    occupy a single common base pair

    :param paths: Analysis path instance
    :param seed: Random number generator seed
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

    # TODO: Finish unit testing


def _non_randomized(paths: AnalysisPaths):
    """
    Implements:

    sort -u $BED > $TMP/$NAME.epitopes.ori2

    :param paths: AnalysisPaths instance
    """
    subprocess_pipes.pipe(
        ["sort", "-u", paths.bed_path], output_path=paths.exclusions_shuffled
    )

    # TODO: Shouldn't this be instead "sort -k 1,1 -k2,2n -u $BED > ..."
    #         sort -k 1,1 -k2,2n -u $BED > ...
    #         Consider the following case:
    #         ENST00000001008 189     198
    #         ENST00000001008 27      36
    #         If sort -u then this is already sorted!
    #         (presuming hierarchy should be chrom, then start, then stop)


class NonRandom(_Randomize):
    """No randomization implemented"""

    @staticmethod
    def apply(params: Parameters):
        _Randomize.check_ready(params)
        _non_randomized(params)


class RandomizeWithoutRegions(_Randomize):
    """Randomizes without user input file"""

    @staticmethod
    def apply(params: Parameters):
        _Randomize.check_ready(params)
        _define_excluded_regions_for_randomization(params)
        _sort_excluded_regions_for_randomization(params, seed=params.seed)


class RandomizeWithRegions(_Randomize):
    """Randomizes with user input file"""

    @staticmethod
    def apply(params: Parameters):
        _Randomize.check_ready(params)
        _randomize_with_target_file(
            params, params.transcripts, seed=params.seed
        )


def _exclude_positively_selected_genes_disabled(paths: AnalysisPaths):
    """
    Implements:

    cp $TMP/$NAME.epitopes.ori2 $TMP/$NAME.epitopes.bed

    No driver genes excluded, so simply copies epitopes file.

    :param paths: AnalysisPaths instance
    """
    subprocess_pipes.pipe(["cp", paths.exclusions_shuffled, paths.epitopes])


def _exclude_positively_selected_genes(
    paths: AnalysisPaths, aux_paths: AuxiliaryPaths
):
    """
    Implements:

    fgrep -w -v -f $SUPA/genes2exclude.txt $TMP/$NAME.epitopes.ori2 >
        $TMP/$NAME.epitopes.bed

    Matches transcript ids not found in genes2exclude.txt to create epitope
    bed file

    :param paths: AnalysisPaths instance
    :param aux_paths: AuxiliaryPaths instance - contains definitions of
            which genes to exclude.
    """

    subprocess_pipes.pipe(
        [
            "fgrep",
            "-w",
            "-v",
            "-f",
            aux_paths.genes_to_exclude.as_posix(),
            paths.exclusions_shuffled.as_posix(),
        ],
        output_path=paths.epitopes,
    )


class _GeneExclusions(_PipelineComponent):
    """Intermediate class for gene exclusions"""

    @staticmethod
    def check_ready(params: Parameters):
        if not params.exclusions_shuffled.exists():
            raise MissingDataError(f"{params.exclusions_shuffled}")


class GeneExclusions(_GeneExclusions):
    """Applies gene exclusions"""

    @staticmethod
    def apply(params: Parameters):
        _GeneExclusions.check_ready(params)
        _exclude_positively_selected_genes(params, AuxiliaryFiles)


class GeneExclusionsDisabled(_GeneExclusions):
    """No gene exclusions"""

    @staticmethod
    def apply(params: Parameters):
        _GeneExclusions.check_ready(params)
        _exclude_positively_selected_genes_disabled(params)


def _get_protein_complement(paths: AnalysisPaths):
    """
    Implements:

    sortBed -i $TMP/$NAME.epitopes.bed |
        complementBed -i stdin -g $TMP/$NAME.protein_length_filt.txt >
            $TMP/$NAME.intra_epitopes_prot.tmp

    cut -f1 $TMP/$NAME.epitopes.bed |
        sort -u |
            fgrep -w -f - $TMP/$NAME.intra_epitopes_prot.tmp >
                $TMP/$NAME.intra_epitopes_prot.bed

    Finds complement to epitope file, subject to filtered protein transcript

    :param paths: AnalysisPaths instance
    """

    subprocess_pipes.pipe(
        ["sortBed", "-i", paths.epitopes],
        [
            "complementBed",
            "-i",
            "stdin",
            "-g",
            paths.filtered_protein_transcript.as_posix(),
        ],
        output_path=paths.intra_epitopes_prot_tmp,
    )

    subprocess_pipes.pipe(
        ["cut", "-f1", paths.epitopes],
        ["sort", "-u"],
        [
            "fgrep",
            "-w",
            "-f",
            "-",
            paths.intra_epitopes_prot_tmp.as_posix(),
        ],
        output_path=paths.intra_epitopes_prot,
    )


class BuildProteinComplement(_PipelineComponent):
    """Build protein complement file"""

    @staticmethod
    def check_ready(params: Parameters):
        for path in (params.epitopes, params.filtered_protein_transcript):
            if not path.exists():
                raise MissingDataError(path)

    @staticmethod
    def apply(params: Parameters):
        BuildProteinComplement.check_ready(params)
        _get_protein_complement(params)


def _prep_ssb192(paths: AnalysisPaths):
    """
    Implements:

    awk '{OFS="\t"}{if( (($2*3)-6) >= 0 )
        {print $1,($2*3)-6,$3*3+3,$0}
            else{print $1,($2*3)-3,$3*3+3,$0}}' $TMP/$NAME.epitopes.bed >
                $TMP/$NAME.epitopes_cds.bed

    :param paths: AnalysisPaths instance
    """

    subprocess_pipes.pipe(
        [
            "awk",
            '{OFS="\t"}{'
            "if( (($2*3)-6) >= 0 )"
            "{print $1,($2*3)-6,$3*3+3,$0}"
            "else{print $1,($2*3)-3,$3*3+3,$0}"
            "}",
            paths.epitopes.as_posix(),
        ],
        output_path=paths.epitopes_cds,
    )


def _prep_not_ssb192(paths: AnalysisPaths):
    """
    Implements:

    awk '{OFS="\t"}{print $1,($2*3)-3,$3*3,$0}' $TMP/$NAME.epitopes.bed >
        $TMP/$NAME.epitopes_cds.bed

    :param paths: AnalysisPaths instance
    """

    subprocess_pipes.pipe(
        ["awk", '{OFS="\t"}{print $1,($2*3)-3,$3*3,$0}', paths.epitopes],
        output_path=paths.epitopes_cds,
    )


class _SSB192Selection(_PipelineComponent):
    """Intermediate class for ssb192 mutrate selection"""

    @staticmethod
    def check_ready(params: Parameters):
        if not params.epitopes.exists():
            raise MissingDataError(params.epitopes.as_posix())


class UseSSB192(_SSB192Selection):
    """Applies ssb192 selection in CDS coordinate prep"""

    @staticmethod
    def apply(params: Parameters):
        _SSB192Selection.check_ready(params)
        _prep_ssb192(params)


class NotSSB192(_SSB192Selection):
    """Does not apply ssb192 selection in CDS coordiante prep"""

    @staticmethod
    def apply(params: Parameters):
        _SSB192Selection.check_ready(params)
        _prep_not_ssb192(params)


def transform_coordinates(paths: AnalysisPaths):
    """
    Implements:

    sortBed -i $TMP/$NAME.epitopes_cds.bed |
        complementBed -i stdin -g $TMP/$NAME.transcript_length_filt.txt >
            $TMP/$NAME.intra_epitopes.tmp

    cut -f1 $TMP/$NAME.epitopes.bed |
        sort -u |
            fgrep -w -f - $TMP/$NAME.intra_epitopes.tmp >
                $TMP/$NAME.intra_epitopes_cds.bed

    Computes the complement for the epitope file in CDS coords.

    :param paths: AnalysisPaths instance
    """

    subprocess_pipes.pipe(
        ["sortBed", "-i", paths.epitopes_cds.as_posix()],
        [
            "complementBed",
            "-i",
            "stdin",
            "-g",
            paths.filtered_transcript.as_posix(),
        ],
        output_path=paths.intra_epitopes_tmp,
    )

    subprocess_pipes.pipe(
        ["cut", "-f1", paths.epitopes.as_posix()],
        ["sort", "-u"],
        ["fgrep", "-w", "-f", "-", paths.intra_epitopes_tmp.as_posix()],
        output_path=paths.intra_epitopes_cds,
    )


class BuildIntraEpitopesCDS(_PipelineComponent):
    """Builds the complement to the epitope in cds coords"""

    @staticmethod
    def check_ready(params: Parameters):
        for path in (
            params.epitopes,
            params.epitopes_cds,
            params.filtered_transcript,
        ):
            if not path.exists():
                raise MissingDataError(path)

    @staticmethod
    def apply(params: Parameters):
        BuildIntraEpitopesCDS.check_ready(params)
        transform_coordinates(params)
