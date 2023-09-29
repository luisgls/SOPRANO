# Installation

SOPRANO requires Mac or Linux OS. All dependencies (including
bioinformatics tools such as VEP and bedtools) are built within a conda
environment for coherent reproducibility. The suite of SOPRANO tools is
assembled as a Python package, which is installed into the encompassing conda
environment.

To install the environment and SOPRANO code, users are therefore required
to have conda (or mamba) available on their system.

**From the repository root, run the command**

```shell
. setup.sh
```

This will automatically construct and activate the `soprano-dev` environment.
After the installation, users can activate the pre-built environment with

```shell
conda activate soprano-dev
```

### Developer notes

If wish to contribute to the development of SOPRANO, run the
setup script with the additional argument

```shell
. setup.sh dev
```

This will install additional packages (such as black, ruff, isort and pytest).
From the repository root, calling `pre-commit install` will prepare a number of
hooks that will enforce code style consistent with the wider repository.
