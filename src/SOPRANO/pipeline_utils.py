import pathlib
from datetime import datetime

from SOPRANO.objects import Parameters
from SOPRANO.prepare_coordinates import (
    _non_randomized,
    filter_transcript_files,
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
    """No randomization implemented"""

    def apply(self, params: Parameters):
        _non_randomized(params)
