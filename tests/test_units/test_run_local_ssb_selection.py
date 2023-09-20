import pathlib

import pytest

from SOPRANO import run_local_ssb_selection as local_ssb


def test_check_path():
    # Path is None checks
    local_ssb.check_path(None, optional=True)
    with pytest.raises(Exception):
        local_ssb.check_path(None, optional=False)

    # Path is not None checks
    local_ssb.check_path(pathlib.Path.cwd())
    with pytest.raises(Exception):
        local_ssb.check_path(pathlib.Path.cwd() / "foo")
