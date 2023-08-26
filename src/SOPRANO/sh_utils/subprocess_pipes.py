import pathlib
import subprocess


def process_output_to_string(ps: subprocess.CompletedProcess) -> str:
    return ps.stdout.decode("utf-8").strip()


def process_output_to_file(
    ps: subprocess.CompletedProcess, path: pathlib.Path, overwrite=False
) -> None:
    """

    :param ps: complete subprocess
    :param path: Path to write stdout to
    :param overwrite: if False, exception raised if output_path exists
    """
    if path.exists() and not overwrite:
        raise FileExistsError(path)

    with open(path, "w") as f:
        f.writelines(process_output_to_string(ps))


def pipe(
    *args: list | tuple, _input=None, output_path=None, overwrite=False
) -> str:
    """
    Recursive pipe of cmdline args through subprocesses
    :param args: comma seperated lists or tuples representing cmds
    :param _input: input to pipe into leading command from args
    :param output_path: Path to write stdout to
    :param overwrite: if False, exception raised if output_path exists
    :return: string representation of stdout from cumulative piped processes
    """

    ps = subprocess.run(args[0], input=_input, capture_output=True)

    if len(args) > 1:
        return pipe(
            *args[1:],
            _input=ps.stdout,
            output_path=output_path,
            overwrite=overwrite,
        )
    else:
        if output_path is not None:
            process_output_to_file(ps, path=output_path, overwrite=overwrite)

        return process_output_to_string(ps)
