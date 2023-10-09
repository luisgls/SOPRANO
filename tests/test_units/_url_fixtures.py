import pickle as pk

import pytest

from SOPRANO.utils.url_utils import _SOPRANO_ENSEMBL_RELEASES


@pytest.fixture
def backup_releases_cache():
    if not _SOPRANO_ENSEMBL_RELEASES.exists():
        with open(_SOPRANO_ENSEMBL_RELEASES, "wb") as f:
            pk.dump({}, f)

    backup_name = _SOPRANO_ENSEMBL_RELEASES.with_suffix(".releases.backup")
    backup_releases_path = _SOPRANO_ENSEMBL_RELEASES.parent.joinpath(
        backup_name
    )
    _SOPRANO_ENSEMBL_RELEASES.rename(backup_releases_path)
    return backup_releases_path
