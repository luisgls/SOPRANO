import pytest

from SOPRANO.env_utils import config_conda_env


@pytest.mark.dependency(name="has_conda")
def test__has_conda():
    assert config_conda_env._has_conda() or config_conda_env._has_mamba()


@pytest.mark.dependency(depends=["has_conda"])
def test_running_soprano_env():
    assert config_conda_env.running_soprano_env()


@pytest.mark.dependency(depends=["has_conda"])
def test_vep_installed():
    assert config_conda_env.vep_installed()


@pytest.mark.dependency(depends=["has_conda"])
def test_bedtools_installed():
    assert config_conda_env.bedtools_installed()
