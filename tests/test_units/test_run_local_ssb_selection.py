import pathlib

import pytest

import SOPRANO.utils.misc_utils


def test_check_path():
    # Path is None checks
    SOPRANO.utils.misc_utils.check_cli_path(None, optional=True)
    with pytest.raises(Exception):
        SOPRANO.utils.misc_utils.check_cli_path(None, optional=False)

    # Path is not None checks
    SOPRANO.utils.misc_utils.check_cli_path(pathlib.Path.cwd())
    with pytest.raises(Exception):
        SOPRANO.utils.misc_utils.check_cli_path(pathlib.Path.cwd() / "foo")
