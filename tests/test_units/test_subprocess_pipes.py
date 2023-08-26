import pathlib
import subprocess
import tempfile

import pytest

from SOPRANO.sh_utils import subprocess_pipes as sbp


def test_process_output_to_string():
    test_process = subprocess.run(["echo", "123"], capture_output=True)

    assert sbp.process_output_to_string(test_process) == "123"


def test_process_output_to_file():
    test_process = subprocess.run(["echo", "123"], capture_output=True)

    with tempfile.TemporaryDirectory() as tmpdir:
        output_path = pathlib.Path(tmpdir).joinpath("out")

        # write subprocess output to file
        sbp.process_output_to_file(test_process, output_path)

        assert output_path.exists()

        # Check error raised when attempt to overwrite
        with pytest.raises(FileExistsError):
            sbp.process_output_to_file(test_process, output_path)

        # Check on error raised when overwrite=True
        sbp.process_output_to_file(test_process, output_path, overwrite=True)

        # Read output as written to file
        with open(output_path, "r") as f:
            lines = f.readlines()

        # Check expectation
        assert len(lines) == 1
        assert lines[0] == "123"
