import pathlib
import subprocess
import tempfile

import pytest

from SOPRANO.utils import sh_utils as sbp


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


def test_pipe():
    with tempfile.TemporaryDirectory() as _tmpdir:
        tmpdir = pathlib.Path(_tmpdir)

        # write files with suffix 00, 05, 10, 15 (for grep tests)
        for idx in range(4):
            file_name = tmpdir.joinpath("file_%02d" % (idx * 5))
            with open(file_name, "w"):
                pass
            assert pathlib.Path(file_name).exists()

        assert sbp.pipe(["ls", "-1", tmpdir], ["grep", "05"]) == "file_05"

        output_path = tmpdir.joinpath("out")

        sbp.pipe(
            ["ls", "-1", tmpdir],
            ["grep", "file"],
            ["grep", "0"],
            output_path=output_path,
        )

        expected_lines = ["file_00", "file_05"]

        with open(output_path, "r") as f:
            written_lines = f.readlines()

        for e, w in zip(expected_lines, written_lines):
            assert e == w.strip()
