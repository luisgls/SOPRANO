import os
import pathlib
import subprocess

SH_UTILS_DIR = pathlib.Path(__file__).resolve().parent
SOPRANO_SRC_DIR = SH_UTILS_DIR.parent

# This should be pointed to within the ~/.condarc file
# TODO: no longer needed...
CONDA_DIR = SOPRANO_SRC_DIR / "conda_env"


def _create_new_condarc(conda_env_location: str, condarc_path: str) -> None:
    # TODO: no longer needed...
    condarc_lines = [
        "envs_dirs:\n",
        "  - ~/.conda/envs\n",
        f"  - {conda_env_location}\n",
    ]

    with open(condarc_path, "w") as f:
        f.writelines(condarc_lines)


def _update_condarc(conda_env_location: str, condarc_path: str) -> None:
    # TODO: no longer needed...
    local_env_in_rc = False
    envs_clause_in_rc = False
    envs_clause = "envs_dirs:\n"
    envs_clause_idx = -1

    # Read lines from condarc file
    with open(condarc_path, "r") as f:
        condarc_lines = f.readlines()

    # Iterate over lines (and count their index)
    for line_idx, line in enumerate(condarc_lines):
        if conda_env_location in line:
            local_env_in_rc = True
        if envs_clause in line:
            envs_clause_in_rc = True
            envs_clause_idx = line_idx

    if not local_env_in_rc:
        if envs_clause_in_rc:
            # Get yml tab spacing from following line
            # crate new line with custom location
            # and insert into list of original lines
            yml_tab_spacing = condarc_lines[envs_clause_idx + 1].split("-")[0]
            new_line = f"{yml_tab_spacing}- {conda_env_location}\n"
            condarc_lines.insert(envs_clause_idx + 2, new_line)
        else:
            # Insert envs clause into start of rc file
            # then insert standard and custom locations line
            std_line = "  - ~/.conda/envs\n"
            new_line = f"  - {conda_env_location}\n"
            condarc_lines.insert(0, envs_clause)
            condarc_lines.insert(1, std_line)
            condarc_lines.insert(2, new_line)

        with open(condarc_path, "w") as f:
            f.writelines(condarc_lines)


def _has_conda():
    which_conda = subprocess.run(["which", "conda"], stdout=subprocess.PIPE)
    return which_conda.returncode == 0


def _has_mamba():
    which_mamba = subprocess.run(["which", "mamba"], stdout=subprocess.PIPE)
    return which_mamba.returncode == 0


def _conda_env_exists(_flavor="conda") -> bool:
    # TODO: no longer needed...
    grep_env = subprocess.run([_flavor, "env", "list"], stdout=subprocess.PIPE)
    return grep_env.returncode == 0


def _update_conda_env(conda_env_location: str):
    # TODO: no longer needed...
    pass


class CondaEnvironmentFailure(Exception):
    pass


def _build_conda_env(conda_env_location: str):
    # TODO: no longer needed...
    yml_path = SOPRANO_SRC_DIR / "local.yml"

    if _conda_env_exists():
        #     _update_conda_env(conda_env_location)
        # else:
        conda_exec = "mamba" if _has_mamba() else "conda"

        soprano_env_location = os.path.join(conda_env_location, "SOPRANO")

        if os.path.exists(soprano_env_location):
            raise OSError(
                f"'SOPRANO' environment directory already exists: "
                f"{soprano_env_location}"
            )

        print("-- Building conda environment, this will take some time...")
        os.environ["CONDA_ALWAYS_YES"] = "true"
        env_build_output = subprocess.run(
            [
                f"{conda_exec}",
                "env",
                "create",
                "--prefix",
                f"{soprano_env_location}",
                "--file",
                f"{yml_path}",
            ],
            stdout=subprocess.PIPE,
        )

        if env_build_output.returncode != 0:
            raise CondaEnvironmentFailure(
                "Non-zero exit code generated from conda environment build."
            )


def prepare_condarc():
    # TODO: no longer needed...
    conda_rc_path = pathlib.Path.home() / ".condarc"

    if not conda_rc_path.exists():
        _create_new_condarc(CONDA_DIR.as_posix(), conda_rc_path.as_posix())
    else:
        _update_condarc(CONDA_DIR.as_posix(), conda_rc_path.as_posix())


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
