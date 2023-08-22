import pathlib

SH_UTILS_DIR = pathlib.Path(__file__).resolve().parent
SOPRANO_SRC_DIR = SH_UTILS_DIR.parent
CONDA_DIR = SOPRANO_SRC_DIR / "conda_env"
