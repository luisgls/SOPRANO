from pathlib import PosixPath

import pytest

from SOPRANO import run_local_ssb_selection as local_ssb


def test__validate_input_path():
    # Path is None checks
    local_ssb._validate_input_path(None, optional=True)
    with pytest.raises(Exception):
        local_ssb._validate_input_path(None, optional=False)

    # Path is not None checks
    local_ssb._validate_input_path(PosixPath.cwd())
    with pytest.raises(Exception):
        local_ssb._validate_input_path(PosixPath.cwd() / "foo")
