import os
import pathlib

# Common dirs from source root src/SOPRANO
_SOPRANO_SRC = pathlib.Path(__file__).parent.parent
_SOPRANO_SCRIPTS = _SOPRANO_SRC / "scripts"
_SOPRANO_R = _SOPRANO_SRC / "R"
_SOPRANO_DATA = _SOPRANO_SRC / "data"
_SOPRANO_IMMUNO = _SOPRANO_SRC / "immunopeptidomes"
_SOPRANO_IMMUNO_HUMANS = _SOPRANO_IMMUNO / "human"
_SOPRANO_EXAMPLES = _SOPRANO_SRC / "examples"
_SOPRANO_INSTALLERS = _SOPRANO_SRC / "shell_utils"

# Common dirs from repository root
_SOPRANO_REPO = _SOPRANO_SRC.parent.parent
_SOPRANO_DEFAULT_CACHE = _SOPRANO_REPO / "pipeline_cache"
_SOPRANO_ENSEMBL_CACHE = _SOPRANO_REPO / "ensembl_downloads"
_SOPRANO_HOMO_SAPIENS = _SOPRANO_ENSEMBL_CACHE / "homo_sapiens"

# Test dirs
_SOPRANO_TESTS = _SOPRANO_REPO / "tests"
_SOPRANO_UNIT_TESTS = _SOPRANO_TESTS / "units"
_SOPRANO_INT_TESTS = _SOPRANO_TESTS / "e2e"
_SOPRANO_CFG_TESTS = _SOPRANO_TESTS / "test_configuration"

# Common dirs from app sources
_SOPRANO_APP_SOURCES = _SOPRANO_REPO / "app_sources"
_SOPRANO_APP_ANNOTATED_INPUTS = _SOPRANO_APP_SOURCES / "annotated_inputs"
_SOPRANO_APP_IMMUNO = _SOPRANO_APP_SOURCES / "immunopeptidomes"
_SOPRANO_APP_COORDS = _SOPRANO_APP_SOURCES / "coordinate_files"

# Other system paths
_STD_SYS_VEP = pathlib.Path.home() / ".vep"


class Directories:
    @staticmethod
    def src(sub_path_item="") -> pathlib.Path:
        return _SOPRANO_SRC.joinpath(sub_path_item)

    @staticmethod
    def scripts(sub_path_item="") -> pathlib.Path:
        return _SOPRANO_SCRIPTS.joinpath(sub_path_item)

    @staticmethod
    def r_scripts(sub_path_item="") -> pathlib.Path:
        return _SOPRANO_R.joinpath(sub_path_item)

    @staticmethod
    def data(sub_path_item="") -> pathlib.Path:
        return _SOPRANO_DATA.joinpath(sub_path_item)

    @staticmethod
    def ensembl_downloads(sub_path_item="") -> pathlib.Path:
        return _SOPRANO_ENSEMBL_CACHE.joinpath(sub_path_item)

    @staticmethod
    def genomes_homo_sapiens(sub_path_item="") -> pathlib.Path:
        return _SOPRANO_HOMO_SAPIENS.joinpath(sub_path_item)

    @staticmethod
    def immunopeptidomes(sub_path_item="") -> pathlib.Path:
        return _SOPRANO_IMMUNO.joinpath(sub_path_item)

    @staticmethod
    def immunopeptidomes_humans(sub_path_item="") -> pathlib.Path:
        return _SOPRANO_IMMUNO_HUMANS.joinpath(sub_path_item)

    @staticmethod
    def examples(sub_path_item="") -> pathlib.Path:
        return _SOPRANO_EXAMPLES.joinpath(sub_path_item)

    @staticmethod
    def cache(sub_path_item="") -> pathlib.Path:
        if "SOPRANO_CACHE" in os.environ.keys():
            active_cache = pathlib.Path(os.environ["SOPRANO_CACHE"])
        else:
            if not _SOPRANO_DEFAULT_CACHE.exists():
                _SOPRANO_DEFAULT_CACHE.mkdir()

            active_cache = _SOPRANO_DEFAULT_CACHE

        return active_cache.joinpath(sub_path_item)

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

    @staticmethod
    def std_sys_vep(sub_path_item="") -> pathlib.Path:
        return _STD_SYS_VEP.joinpath(sub_path_item)

    @staticmethod
    def app_sources(sub_path_item="") -> pathlib.Path:
        return _SOPRANO_APP_SOURCES.joinpath(sub_path_item)

    @staticmethod
    def app_annotated_inputs(sub_path_item="") -> pathlib.Path:
        return _SOPRANO_APP_ANNOTATED_INPUTS.joinpath(sub_path_item)

    @staticmethod
    def app_coordinate_files(sub_path_item="") -> pathlib.Path:
        return _SOPRANO_APP_COORDS.joinpath(sub_path_item)

    @staticmethod
    def app_immunopeptidomes(sub_path_item="") -> pathlib.Path:
        return _SOPRANO_APP_IMMUNO.joinpath(sub_path_item)


def is_empty(path: pathlib.Path) -> bool:
    """
    Checks whether file at path has size of zero
    :param path: pathlib Path object
    :return: True if path is empty else False
    """
    return path.stat().st_size == 0


def _check_paths(*dependent_paths: pathlib.Path):
    for path in dependent_paths:
        if not path.exists():
            raise FileNotFoundError(path)


def check_cli_path(cli_path: pathlib.Path | None, optional=False):
    if cli_path is None:
        if not optional:
            raise FileNotFoundError(
                "Input path is not optional and path is None!"
            )
    elif not cli_path.exists():
        raise FileNotFoundError(f"CLI input path does not exist: {cli_path}")


def genome_pars_to_paths(ref, release):
    """
    Translates human genome reference and release ids into a tuple of paths

    - genome reference fasta file path
    - chrom sizes file path

    :param ref: Genome reference ID
    :param release: Ensembl release ID
    :return: Tuple of paths: reference fasta file, chrom sizes
    """
    data_dir = Directories.genomes_homo_sapiens(f"{release}_{ref}")
    genome_path = data_dir.joinpath(f"Homo_sapiens.{ref}.dna.toplevel.fa")
    chroms_path = data_dir.joinpath(f"Homo_sapiens.{ref}.dna.toplevel.chrom")
    return genome_path, chroms_path
