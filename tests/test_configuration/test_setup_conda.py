import pytest

from SOPRANO.utils.env_utils import (
    _has_conda,
    _has_mamba,
    bedtools_installed,
    running_soprano_env,
    vep_installed,
)


@pytest.mark.dependency(name="has_conda")
def test__has_conda():
    assert _has_conda() or _has_mamba()


@pytest.mark.dependency(depends=["has_conda"])
def test_running_soprano_env():
    assert running_soprano_env()


@pytest.mark.dependency(depends=["has_conda"])
def test_vep_installed():
    assert vep_installed()


@pytest.mark.dependency(depends=["has_conda"])
def test_bedtools_installed():
    assert bedtools_installed()
