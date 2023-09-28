import os
import pathlib

_SOPRANO_SRC = pathlib.Path(__file__).parent
_SOPRANO_SCRIPTS = _SOPRANO_SRC / "scripts"
_SOPRANO_DATA = _SOPRANO_SRC / "data"
_SOPRANO_HOMO_SAPIENS = _SOPRANO_DATA / "homo_sapiens"
_SOPRANO_IMMUNO = _SOPRANO_SRC / "immunopeptidomes"
_SOPRANO_HUMANS = _SOPRANO_IMMUNO / "human"
_SOPRANO_EXAMPLES = _SOPRANO_SRC / "examples"
_SOPRANO_REPO = _SOPRANO_SRC.parent.parent
_SOPRANO_DEFAULT_CACHE = _SOPRANO_REPO / "pipeline_cache"
_SOPRANO_TESTS = _SOPRANO_REPO / "tests"
_SOPRANO_UNIT_TESTS = _SOPRANO_TESTS / "test_units"
_SOPRANO_INT_TESTS = _SOPRANO_TESTS / "test_integrations"
_SOPRANO_CFG_TESTS = _SOPRANO_TESTS / "test_configuration"
_SOPRANO_INSTALLERS = _SOPRANO_SRC / "bash_installers"


class Directories:
    @staticmethod
    def src(sub_path_item="") -> pathlib.Path:
        return _SOPRANO_SRC.joinpath(sub_path_item)

    @staticmethod
    def scripts(sub_path_item="") -> pathlib.Path:
        return _SOPRANO_SCRIPTS.joinpath(sub_path_item)

    @staticmethod
    def data(sub_path_item="") -> pathlib.Path:
        return _SOPRANO_DATA.joinpath(sub_path_item)

    @staticmethod
    def homo_sapien_genomes(sub_path_item="") -> pathlib.Path:
        return _SOPRANO_HOMO_SAPIENS.joinpath(sub_path_item)

    @staticmethod
    def immunopeptidomes(sub_path_item="") -> pathlib.Path:
        return _SOPRANO_IMMUNO.joinpath(sub_path_item)

    @staticmethod
    def immuno_humans(sub_path_item="") -> pathlib.Path:
        return _SOPRANO_HUMANS.joinpath(sub_path_item)

    @staticmethod
    def examples(sub_path_item="") -> pathlib.Path:
        return _SOPRANO_EXAMPLES.joinpath(sub_path_item)

    @staticmethod
    def cache(sub_path_item="") -> pathlib.Path:
        if "SOPRANO_CACHE" in os.environ.keys():
            return pathlib.Path(os.environ["SOPRANO_CACHE"])

        if not _SOPRANO_DEFAULT_CACHE.exists():
            _SOPRANO_DEFAULT_CACHE.joinpath(sub_path_item).mkdir()

        return _SOPRANO_DEFAULT_CACHE.joinpath(sub_path_item)

    @staticmethod
    def tests(sub_path_item="") -> pathlib.Path:
        return _SOPRANO_TESTS.joinpath(sub_path_item)

    @staticmethod
    def unit_tests(sub_path_item="") -> pathlib.Path:
        return _SOPRANO_UNIT_TESTS.joinpath(sub_path_item)

    @staticmethod
    def int_tests(sub_path_item="") -> pathlib.Path:
        return _SOPRANO_INT_TESTS.joinpath(sub_path_item)

    @staticmethod
    def cfg_tests(sub_path_item="") -> pathlib.Path:
        return _SOPRANO_CFG_TESTS.joinpath(sub_path_item)

    @staticmethod
    def installers(sub_path_item="") -> pathlib.Path:
        return _SOPRANO_INSTALLERS.joinpath(sub_path_item)


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
