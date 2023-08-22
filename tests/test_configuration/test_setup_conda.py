import pathlib
import tempfile

from SOPRANO.sh_utils import setup_conda


def test_running_conda():
    """
    TODO: Check whether python interpreter is linked to conda env
    :return:
    """
    pass


def test__create_new_condarc():
    with tempfile.TemporaryDirectory() as tmp_dir:
        condarc_path = pathlib.Path(tmp_dir).joinpath(".condarc").as_posix()
        mock_env_location = pathlib.Path.cwd().as_posix()

        expected_lines = [
            "envs_dirs:\n",
            "  - ~/.conda/envs\n",
            f"  - {mock_env_location}\n",
        ]

        setup_conda._create_new_condarc(mock_env_location, condarc_path)
        with open(condarc_path, "r") as f:
            written_lines = f.readlines()

        for expected, written in zip(expected_lines, written_lines):
            assert expected == written, (expected, written)


def test__update_condarc():
    assert False


def test_prepare_condarc():
    assert False
