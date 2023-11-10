# Installation

### Dependencies

Though SOPRANO has been developed as a Python package, it relies on a variety
of languages and tools:

- `Python >= 3.10`
- `R >= 4.0.0`
- `Perl 5`
- `Bedtools >= 5.32.0`

### Conda environment

The most straightforward way to manage these dependencies is to create a
`conda` environment. From the repository root:

```shell
conda env create -f env.yml
```

After creation, the environment can be activated with

```shell
conda activate soprano
```

### Non-conda users

Users who do not wish to use `conda` should ensure that they have appropriate
versions of `Python`, `R`, `Perl` and `Bedtools` installed on their system.

To install SOPRANO as a Python package, it is recommended to create a Python
virtual environment, e.g., from the repository root

```shell
python3 -m venv .venv
```

After creation, the environment can be activated with

```shell
source .venv/bin/activate
```

### Installing R package dependencies

To install the dependencies for mutation annotations, from the repository
root

```shell
Rscript install_R_pkgs.R
```

Conda users should readily have most of these packages installed from the
previous environment build, though `dndscv` will still be installed directly
from GitHub.

### Installing the SOPRANO Python package

Ensure an appropriate environment is activated, then

```shell
pip install -e .
```

This will install a standard (though editable version of SOPRANO).

### Validating the installation

To check that the validity of your installation, users can run

```shell
pytest -s tests/test_installation.py
```

### Developer notes

Developers should install the SOPRANO package with additional dependencies

```shell
pip install -e .[dev]
```

This will install additional packages (such as black, ruff, isort and pytest).
Consistent styling can be enforced with pre-commit hooks: `pre-commit install`.
