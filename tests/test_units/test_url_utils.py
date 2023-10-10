import pickle as pk

import pytest

import SOPRANO.core.objects
from SOPRANO.utils import url_utils


def test_filename_from_url():
    base = "path/to/dir"
    file_name = "something.gz"

    joined = base + "/" + file_name

    assert url_utils.filename_from_url(joined) == file_name


def test_download_from_url(tmp_path):
    test_url = "https://www.icr.ac.uk/assets/img/logo_share.jpg"
    test_target_path = tmp_path.joinpath("icr_logo.jpg")
    url_utils.download_from_url(test_url, target_path=test_target_path)
    assert test_target_path.exists()


def test_get_dna_url():
    species = "foo"
    dna_url = url_utils.get_dna_url(species)
    assert dna_url.endswith(f"{species}/dna")


def test_check_ensembl_dna_url():
    dna_url = (
        "https://ftp.ensembl.org/pub/release-{RELEASE}/"
        "fasta/homo_sapiens/dna"
    )

    url_utils.check_ensembl_dna_url(dna_url, 110)  # No error expected

    with pytest.raises(url_utils.BadResponse):
        url_utils.check_ensembl_dna_url(dna_url, 1000)  # Error expected


def test_check_ensembl_file_url():
    toplevel_url = (
        "https://ftp.ensembl.org/pub/grch37/release-{RELEASE}/"
        "fasta/homo_sapiens/dna/"
        "Homo_sapiens.GRCh37.dna.toplevel.fa.gz"
    )

    url_utils.check_ensembl_file_url(toplevel_url, 110)  # No error expected

    with pytest.raises(url_utils.BadResponse):
        url_utils.check_ensembl_file_url(toplevel_url, 1000)


def test_get_cached_release_dict():
    retrieved = url_utils.get_cached_release_dict()
    assert isinstance(retrieved, dict)


def test_get_cached_release_value(backup_releases_cache):
    backup_releases_path = backup_releases_cache

    url_to_return_defaults = "https://pytest4defaults.soprano"
    assert [-1, -1] == url_utils.get_cached_release_value(
        url_to_return_defaults
    )

    url_to_return_values = "https://pytest4values.soprano"
    url_values = [100, 105]
    url_dict = {url_to_return_values: url_values}

    with open(url_utils._SOPRANO_ENSEMBL_RELEASES, "wb") as f:
        pk.dump(url_dict, f)

    assert url_values == url_utils.get_cached_release_value(
        url_to_return_values
    )

    url_utils._SOPRANO_ENSEMBL_RELEASES.unlink(missing_ok=False)
    backup_releases_path.rename(url_utils._SOPRANO_ENSEMBL_RELEASES)


def test_find_earliest_release(backup_releases_cache):
    toplevel_url = (
        "https://ftp.ensembl.org/pub/release-{RELEASE}/fasta/"
        "homo_sapiens/dna/Homo_sapiens.GRCh38.dna.toplevel.fa.gz"
    )
    assert url_utils.find_earliest_release(toplevel_url) == 76
    backup_releases_cache.rename(url_utils._SOPRANO_ENSEMBL_RELEASES)


def test_find_latest_release(backup_releases_cache):
    toplevel_url = (
        "https://ftp.ensembl.org/pub/release-{RELEASE}/fasta/"
        "homo_sapiens/dna/Homo_sapiens.GRCh38.dna.toplevel.fa.gz"
    )

    # As of 06 Oct 2023, latest is 110

    # Insert "old" value that is to be updated in cache
    url_utils.cache_earliest_release(toplevel_url, 109)

    # Check that updated value is indeed larger
    assert url_utils.find_latest_release(toplevel_url) > 109
    backup_releases_cache.rename(url_utils._SOPRANO_ENSEMBL_RELEASES)


def test__cache_release_data(backup_releases_cache):
    toplevel_url = "https://ftp.ensembl.org/foo/bar"

    url_utils._cache_release_data(toplevel_url, 100, True)
    assert url_utils.get_cached_release_value(toplevel_url) == [100, -1]

    url_utils._cache_release_data(toplevel_url, 99, True)
    assert url_utils.get_cached_release_value(toplevel_url) == [100, -1]

    url_utils._cache_release_data(toplevel_url, 101, True)
    assert url_utils.get_cached_release_value(toplevel_url) == [101, -1]

    url_utils._cache_release_data(toplevel_url, 110, False)
    assert url_utils.get_cached_release_value(toplevel_url) == [101, 110]

    backup_releases_cache.rename(url_utils._SOPRANO_ENSEMBL_RELEASES)


def test_cache_latest_release(backup_releases_cache):
    toplevel_url = "https://ftp.ensembl.org/foo/bar"

    url_utils._cache_release_data(toplevel_url, 100, False)
    with open(url_utils._SOPRANO_ENSEMBL_RELEASES, "rb") as f:
        _x = pk.load(f)

    url_utils._SOPRANO_ENSEMBL_RELEASES.unlink()

    url_utils.cache_latest_release(toplevel_url, 100)
    with open(url_utils._SOPRANO_ENSEMBL_RELEASES, "rb") as f:
        x = pk.load(f)

    assert _x == x

    backup_releases_cache.rename(url_utils._SOPRANO_ENSEMBL_RELEASES)


def test_cache_earliest_release(backup_releases_cache):
    toplevel_url = "https://ftp.ensembl.org/foo/bar"

    url_utils._cache_release_data(toplevel_url, 100, True)
    with open(url_utils._SOPRANO_ENSEMBL_RELEASES, "rb") as f:
        _x = pk.load(f)

    url_utils._SOPRANO_ENSEMBL_RELEASES.unlink()

    url_utils.cache_earliest_release(toplevel_url, 100)
    with open(url_utils._SOPRANO_ENSEMBL_RELEASES, "rb") as f:
        x = pk.load(f)

    assert _x == x

    backup_releases_cache.rename(url_utils._SOPRANO_ENSEMBL_RELEASES)


def test_build_ensembl_urls():
    species = "foo"
    reference = "bar"

    expected_toplevel = (
        "https://ftp.ensembl.org/pub/release-{RELEASE}/"
        "fasta/foo/dna/"
        "Foo.bar.dna.toplevel.fa.gz"
    )
    expected_primary_assembly = (
        "https://ftp.ensembl.org/pub/release-{RELEASE}/"
        "fasta/foo/dna/"
        "Foo.bar.dna.primary_assembly.fa.gz"
    )

    expected_dict = {
        "toplevel": expected_toplevel,
        "primary_assembly": expected_primary_assembly,
    }

    assert expected_dict == url_utils.build_ensembl_urls(species, reference)


def test__dest_directory(foo_bar_reference):
    x: SOPRANO.core.objects._GatherReferences = foo_bar_reference
    release = 100
    posix_path = x._dest_directory(release).as_posix()

    assert posix_path.split("/")[-2] == x.species
    assert posix_path.split("/")[-1] == f"{release}_{x.reference}"


def test__dest_fa_gz(foo_bar_reference):
    x: SOPRANO.core.objects._GatherReferences = foo_bar_reference
    release = 100

    assert (
        x._dest_fa_gz(release, _toplevel=True).name
        == x.toplevel_url.split("/")[-1]
    )
    assert (
        x._dest_fa_gz(release, _toplevel=False).name
        == x.primary_assembly_url.split("/")[-1]
    )


def test__dest_fa(foo_bar_reference):
    x: SOPRANO.core.objects._GatherReferences = foo_bar_reference
    release = 100

    assert x._dest_fa(release, _toplevel=True).name == x.toplevel_url.split(
        "/"
    )[-1].rstrip(".gz")
    assert x._dest_fa(
        release, _toplevel=False
    ).name == x.primary_assembly_url.split("/")[-1].rstrip(".gz")


def test__dest_chrom(foo_bar_reference):
    x: SOPRANO.core.objects._GatherReferences = foo_bar_reference
    release = 100

    assert x.dest_chrom(release, _toplevel=True).name == x.toplevel_url.split(
        "/"
    )[-1].replace(".fa.gz", ".chrom")


def test_toplevel_fa_gz_path(foo_bar_reference):
    x: SOPRANO.core.objects._GatherReferences = foo_bar_reference
    release = 100

    assert x.toplevel_fa_gz_path(release) == x._dest_fa_gz(
        release, _toplevel=True
    )


def test_toplevel_fa_path(foo_bar_reference):
    x: SOPRANO.core.objects._GatherReferences = foo_bar_reference
    release = 100

    assert x.toplevel_fa_path(release) == x._dest_fa(release, _toplevel=True)


def test_primary_assembly_fa_gz_path(foo_bar_reference):
    x: SOPRANO.core.objects._GatherReferences = foo_bar_reference
    release = 100

    assert x.primary_assembly_fa_gz_path(release) == x._dest_fa_gz(
        release, _toplevel=False
    )


def test_primary_assembly_fa_path(foo_bar_reference):
    x: SOPRANO.core.objects._GatherReferences = foo_bar_reference
    release = 100

    assert x.primary_assembly_fa_path(release) == x._dest_fa(
        release, _toplevel=False
    )


def test__check_release_ok():
    for x in (
        SOPRANO.core.objects.EnsemblData.homo_sapiens_GRCh38(),
        SOPRANO.core.objects.EnsemblData.homo_sapiens_GRCh37(),
    ):
        x._check_release_ok(100)  # Should be fine

        # Check val error raised on out of range releases
        with pytest.raises(ValueError):
            x._check_release_ok(1)
        with pytest.raises(ValueError):
            x._check_release_ok(1000)
