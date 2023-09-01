import pathlib
import subprocess


def process_output_to_string(ps: subprocess.CompletedProcess) -> str:
    return ps.stdout.decode("utf-8").strip()


def process_output_to_file(
    ps: subprocess.CompletedProcess,
    path: pathlib.Path,
    overwrite=False,
    mode="w",
) -> None:
    """

    :param ps: complete subprocess
    :param path: Path to write stdout to
    :param overwrite: if False, exception raised if output_path exists
    :param mode: writing mode (w or a)
    """
    if path.exists() and not overwrite:
        raise FileExistsError(path)

    if mode == "a" and path.exists():
        with open(path, mode="r") as _:
            _lines = _.readlines()
        prepend_newline = not _lines[-1].endswith("\n")
    else:
        prepend_newline = False

    with open(path, mode=mode) as f:
        if prepend_newline:
            f.writelines("\n")
        f.writelines(process_output_to_string(ps))


def pipe(
    *args: list | tuple,
    _input=None,
    output_path=None,
    overwrite=False,
    mode="w",
) -> str:
    """
    Recursive pipe of cmdline args through subprocesses
    :param args: comma seperated lists or tuples representing cmds
    :param _input: input to pipe into leading command from args
    :param output_path: Path to write stdout to
    :param overwrite: if False, exception raised if output_path exists
    :param mode: writing mode for output_path (default is w). Use 'a' to append
    :return: string representation of stdout from cumulative piped processes
    """

    # TODO: Should we be raising an error on non-zero exit status?
    ps = subprocess.run(args[0], input=_input, capture_output=True)

    # TODO: Should we be catching warnings raised by bedtools?
    #       e.g., these have the prefix ***** WARNING
    #       and usually means something stupid has been done
    #       could check stdout for this, and raise BedtoolsWarning?

    if len(args) > 1:
        return pipe(
            *args[1:],
            _input=ps.stdout,
            output_path=output_path,
            overwrite=overwrite,
            mode=mode,
        )
    else:
        if output_path is not None:
            process_output_to_file(
                ps, path=output_path, overwrite=overwrite, mode=mode
            )

        return process_output_to_string(ps)
