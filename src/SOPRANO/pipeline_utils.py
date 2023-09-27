import pathlib
from datetime import datetime

from SOPRANO.objects import Parameters


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

    def check_ready(self, *params: Parameters) -> None:
        pass
