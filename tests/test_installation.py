from SOPRANO.utils.env_utils import (
    app_installed,
    bedtools_installed,
    cli_installed,
    get_genome_installed,
    link_vep_installed,
    perl_installed,
    rscript_installed,
)


def test_bedtools_installed():
    assert bedtools_installed()


def test_rscript_installed():
    assert rscript_installed()


def test_perl_installed():
    assert perl_installed()


def test_app_installed():
    assert app_installed()


def test_cli_installed():
    assert cli_installed()


def test_link_vep_installed():
    assert link_vep_installed()


def test_get_genome_installed():
    assert get_genome_installed()
