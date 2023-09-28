import os
import pathlib
import subprocess


def _has_conda():
    which_conda = subprocess.run(["which", "conda"], stdout=subprocess.PIPE)
    return which_conda.returncode == 0


def _has_mamba():
    which_mamba = subprocess.run(["which", "mamba"], stdout=subprocess.PIPE)
    return which_mamba.returncode == 0


def running_soprano_env():
    try:
        return os.environ["CONDA_DEFAULT_ENV"] == "soprano-dev"
    except KeyError:
        print("-- Warning: No conda env de")
        return False


def _check_exec(exec_name):
    """
    Check executable is found and inside of conda environment
    :param exec_name: name of executable, e.g. vep
    :return: True if found and inside current conda env, else False
    """

    # Initial check on exit code for executable
    which_exec = subprocess.run(["which", exec_name], stdout=subprocess.PIPE)
    if which_exec.returncode != 0:
        return False

    # Build expected path
    conda_env_path = os.environ["CONDA_PREFIX"]
    expected_path = pathlib.Path(conda_env_path).joinpath("bin")
    expected_path = expected_path.joinpath(exec_name)

    # Decode detected path
    detected_path = which_exec.stdout.decode()
    detected_path = detected_path

    return expected_path.as_posix().strip("\n") == detected_path.strip("\n")


def vep_installed():
    return _check_exec("vep")


def bedtools_installed():
    return _check_exec("bedtools")
