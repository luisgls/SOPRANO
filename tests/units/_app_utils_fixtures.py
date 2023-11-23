import pytest

from SOPRANO.utils.path_utils import Directories


@pytest.fixture
def mock_genome_dir():
    assembly, release = "FOO", 999
    mock_dir = (
        Directories.ensembl_downloads("homo_sapiens") / f"{release}_{assembly}"
    )

    return assembly, release, mock_dir
