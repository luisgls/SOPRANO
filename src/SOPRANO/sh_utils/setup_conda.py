import pathlib

SH_UTILS_DIR = pathlib.Path(__file__).resolve().parent
SOPRANO_SRC_DIR = SH_UTILS_DIR.parent
CONDA_DIR = SOPRANO_SRC_DIR / "conda_env"


def _create_new_condarc(conda_env_location: str, condarc_path: str) -> None:
    condarc_lines = [
        "envs_dirs:\n",
        "  - ~/.conda/envs\n",
        f"  - {conda_env_location}\n",
    ]

    with open(condarc_path, "w") as f:
        f.writelines(condarc_lines)


def _update_condarc(conda_env_location: str, condarc_path: str) -> None:
    local_env_in_rc = False
    envs_clause_in_rc = False
    envs_clause = "envs_dirs:"
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
            new_line = f"{yml_tab_spacing}- {conda_env_location}"
            condarc_lines.insert(envs_clause_idx + 2, new_line)
        else:
            # Insert envs clause into start of rc file
            # then insert standard and custom locations line
            std_line = "  - ~/.conda/envs"
            new_line = f"  - {conda_env_location}"
            condarc_lines.insert(0, envs_clause)
            condarc_lines.insert(1, std_line)
            condarc_lines.insert(2, new_line)


def prepare_condarc():
    conda_rc_path = pathlib.Path.home() / ".condarc"

    if not conda_rc_path.exists():
        _create_new_condarc(CONDA_DIR.as_posix(), conda_rc_path.as_posix())
    else:
        _update_condarc(CONDA_DIR.as_posix(), conda_rc_path.as_posix())


# prepare_condarc()
