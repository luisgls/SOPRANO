import pathlib
import subprocess


def process_output_to_string(ps: subprocess.CompletedProcess) -> str:
    return ps.stdout.decode("utf-8").strip()


def process_output_to_file(
    ps: subprocess.CompletedProcess, path: pathlib.Path, overwrite=False
):
    if path.exists() and not overwrite:
        raise FileExistsError(path)

    with open(path, "w") as f:
        f.writelines(process_output_to_string(ps))


def pipe(
    args_1: list | tuple,
    args_2: list | tuple,
    return_output=False,
    output_path=None,
    overwrite=False,
):
    ps_1 = subprocess.run(args_1, capture_output=True)
    ps_2 = subprocess.run(args_2, input=ps_1.stdout, capture_output=True)

    if output_path is not None:
        process_output_to_file(ps_2, path=output_path, overwrite=overwrite)

    if return_output:
        return process_output_to_string(ps_2)


def rec_pipe(
    *args: list | tuple,
    _input=None,
    return_output=False,
    output_path=None,
    overwrite=False,
):
    """
    Recursive pipe of cmdline args through subprocesses
    :param args: comma seperated lists or tuples representing cmds
    :param _input: input to pipe into leading command from args
    :return:
    """

    ps = subprocess.run(args[0], input=_input, capture_output=True)

    print(args[0], ps.stdout)

    if len(args) > 1:
        return rec_pipe(*args[1:], _input=ps.stdout)
    else:
        if output_path is not None:
            process_output_to_file(ps, path=output_path, overwrite=overwrite)

        if return_output:
            return process_output_to_string(ps)


# rec_pipe(["ls", "-l"], ["grep", "xr"], ["grep", "src"])
