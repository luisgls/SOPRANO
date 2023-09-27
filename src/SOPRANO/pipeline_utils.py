import pathlib
from datetime import datetime

from SOPRANO.objects import AuxiliaryFiles, Parameters
from SOPRANO.obtain_fasta_regions import (
    _get_non_target_regions,
    _get_target_fasta_regions,
)
from SOPRANO.prepare_coordinates import (
    _define_excluded_regions_for_randomization,
    _exclude_positively_selected_genes,
    _exclude_positively_selected_genes_disabled,
    _get_protein_complement,
    _non_randomized,
    _prep_not_ssb192,
    _prep_ssb192,
    _randomize_with_target_file,
    _sort_excluded_regions_for_randomization,
    filter_transcript_files,
    transform_coordinates,
)


def time_output():
    now = datetime.now()
    return now.strftime("%d/%m/%Y %H:%M:%S")


def task_output(msg):
    print(f"[{time_output()}] {msg}")


def is_empty(path: pathlib.Path) -> bool:
    """
    Checks whether file at path has size of zero
    :param path: pathlib Path object
    :return: True if path is empty else False
    """
    return path.stat().st_size == 0


class MissingDataError(Exception):
    pass


class SOPRANOError(Exception):
    pass


def _check_paths(*dependent_paths: pathlib.Path):
    for path in dependent_paths:
        if not path.exists():
            raise MissingDataError(path)


class _PipelineComponent:
    """
    Components of the pipeline are designed to follow the pattern:

    Component.apply(params)

    where apply should include the call to check_read()
    to permit execution.

    Pipeline components should override these methods.
    """

    @staticmethod
    def apply(params: Parameters):
        pass

    @staticmethod
    def check_ready(params: Parameters):
        pass


class _PipelineComponent2:
    msg = ""  # Will be printed to stdout with date/time stamp

    def apply(self, params: Parameters):
        task_output(self.msg)
        self.check_ready(params)
        self._apply(params)

    def _apply(self, params: Parameters):
        pass

    def check_ready(self, params: Parameters):
        pass


class FilterTranscripts(_PipelineComponent):
    """
    Filter transcript files with respect to input bed file
    """

    @staticmethod
    def apply(params: Parameters):
        FilterTranscripts.check_ready(params)
        filter_transcript_files(params, params.transcripts)

    @staticmethod
    def check_ready(params: Parameters):
        for path in (
            params.bed_path,
            params.transcripts.transcript_length,
            params.transcripts.protein_transcript_length,
        ):
            if not path.exists():
                raise MissingDataError(path)


class FilterTranscripts2(_PipelineComponent2):
    msg = "Filtering transcripts"

    def _apply(self, params: Parameters):
        filter_transcript_files(params, params.transcripts)

    def check_ready(self, params: Parameters):
        _check_paths(
            params.bed_path,
            params.transcripts.transcript_length,
            params.transcripts.protein_transcript_length,
        )


class _Randomize(_PipelineComponent):
    """Intermediate class for randomization procedures"""

    @staticmethod
    def check_ready(params: Parameters):
        for path in (
            params.filtered_transcript,
            params.filtered_protein_transcript,
        ):
            if not path.exists():
                raise MissingDataError(
                    f"Filtered transcript not found: {path}"
                )


class _Randomize2(_PipelineComponent2):
    """Intermediate class for randomization procedures"""

    def check_ready(self, params: Parameters):
        _check_paths(
            params.filtered_transcript,
            params.filtered_protein_transcript,
        )


class NonRandom(_Randomize):
    """No randomization implemented"""

    @staticmethod
    def apply(params: Parameters):
        _Randomize.check_ready(params)
        _non_randomized(params)


class NonRandom2(_Randomize2):
    msg = "No randomization selected; selecting unique items from bed file"

    def _apply(self, params: Parameters):
        _non_randomized(params)


class RandomizeWithoutRegions(_Randomize):
    """Randomizes without user input file"""

    @staticmethod
    def apply(params: Parameters):
        _Randomize.check_ready(params)
        _define_excluded_regions_for_randomization(params)
        _sort_excluded_regions_for_randomization(params, seed=params.seed)


class RandomizeWithoutRegions2(_Randomize2):
    msg = "Performing randomization without supplement bed file definitions"

    def _apply(self, params: Parameters):
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


class RandomizeWithRegions2(_Randomize2):
    msg = "Performing randomization using supplement bed file definitions"

    def _apply(self, params: Parameters):
        _randomize_with_target_file(
            params, params.transcripts, seed=params.seed
        )


class _GeneExclusions(_PipelineComponent):
    """Intermediate class for gene exclusions"""

    @staticmethod
    def check_ready(params: Parameters):
        if not params.exclusions_shuffled.exists():
            raise MissingDataError(f"{params.exclusions_shuffled}")


class _GeneExclusions2(_PipelineComponent2):
    def check_ready(self, params: Parameters):
        _check_paths(params.exclusions_shuffled)


class GeneExclusions(_GeneExclusions):
    """Applies gene exclusions"""

    @staticmethod
    def apply(params: Parameters):
        _GeneExclusions.check_ready(params)
        _exclude_positively_selected_genes(params, AuxiliaryFiles)


class GeneExclusions2(_GeneExclusions2):
    msg = "Excluding positively selected genes"

    def _apply(self, params: Parameters):
        _exclude_positively_selected_genes(params, AuxiliaryFiles)


class GeneExclusionsDisabled(_GeneExclusions):
    """No gene exclusions"""

    @staticmethod
    def apply(params: Parameters):
        _GeneExclusions.check_ready(params)
        _exclude_positively_selected_genes_disabled(params)


class GeneExclusionsDisabled2(_GeneExclusions2):
    msg = "Retaining positively selected genes"

    def _apply(self, params: Parameters):
        _exclude_positively_selected_genes_disabled(params)


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


class BuildProteinComplement2(_PipelineComponent2):
    msg = "Building protein complement"

    def check_ready(self, params: Parameters):
        _check_paths(params.epitopes, params.filtered_protein_transcript)

    def _apply(self, params: Parameters):
        _get_protein_complement(params)


class _SSB192Selection(_PipelineComponent):
    """Intermediate class for ssb192 mutrate selection"""

    @staticmethod
    def check_ready(params: Parameters):
        if not params.epitopes.exists():
            raise MissingDataError(params.epitopes.as_posix())


class _SSB192Selection2(_PipelineComponent2):
    """Intermediate class for ssb192 mutrate selection"""

    def check_ready(self, params: Parameters):
        _check_paths(params.epitopes)


class UseSSB192(_SSB192Selection):
    """Applies ssb192 selection in CDS coordinate prep"""

    @staticmethod
    def apply(params: Parameters):
        _SSB192Selection.check_ready(params)
        _prep_ssb192(params)


class UseSSB1922(_SSB192Selection2):
    msg = "Preparing CDS coordinates using SSB 192 substitution"

    def _apply(self, params: Parameters):
        _prep_ssb192(params)


class NotSSB192(_SSB192Selection):
    """Does not apply ssb192 selection in CDS coordiante prep"""

    @staticmethod
    def apply(params: Parameters):
        _SSB192Selection.check_ready(params)
        _prep_not_ssb192(params)


class NotSSB1922(_SSB192Selection2):
    msg = "Preparing CDS coordinates using SSB 7 substitution"

    def _apply(self, params: Parameters):
        _prep_not_ssb192(params)


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


class BuildIntraEpitopesCDS2(_PipelineComponent2):
    msg = "Building intra epitope file in CDS coordinates"

    def check_ready(self, params: Parameters):
        _check_paths(
            params.epitopes,
            params.epitopes_cds,
            params.filtered_transcript,
        )

    def _apply(self, params: Parameters):
        transform_coordinates(params)


class ObtainFastaRegions(_PipelineComponent):
    @staticmethod
    def check_ready(params: Parameters):
        paths = (
            params.transcripts.transcript_fasta,
            params.epitopes_cds,
            params.intra_epitopes_cds,
        )

        for path in paths:
            if not path.exists():
                raise MissingDataError(path)

    @staticmethod
    def apply(params: Parameters):
        ObtainFastaRegions.check_ready(params)
        _get_target_fasta_regions(params, params.transcripts)
        _get_non_target_regions(params, params.transcripts)


class ObtainFastaRegions2(_PipelineComponent2):
    msg = "Obtaining fasta regions"

    def check_ready(self, params: Parameters):
        _check_paths(
            params.transcripts.transcript_fasta,
            params.epitopes_cds,
            params.intra_epitopes_cds,
        )

    def _apply(self, params: Parameters):
        _get_target_fasta_regions(params, params.transcripts)
        _get_non_target_regions(params, params.transcripts)
