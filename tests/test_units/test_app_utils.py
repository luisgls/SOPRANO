import os
from pathlib import Path

import pytest

from SOPRANO.core.objects import GenomePaths
from SOPRANO.utils.app_utils import (
    PipelineUIOptions,
    PipelineUIProcessing,
    _select_from_dict,
)
from SOPRANO.utils.path_utils import _SOPRANO_DEFAULT_CACHE, Directories


class TemporaryFiles:
    def __init__(self, pass_ext, fail_ext, target_dir, mkdir, files):
        self._pass_ext = pass_ext
        self._fail_ext = fail_ext
        self._target_dir = target_dir
        self._mkdir = mkdir
        self._files = files

    @classmethod
    def tmp_in_dir(cls, target_dir: Path, ext: str):
        files = target_dir / f"pytest.{ext}", target_dir / "pytest.fail"
        return cls(
            pass_ext=files[0],
            fail_ext=files[1],
            target_dir=target_dir,
            mkdir=False,
            files=files,
        )

    @classmethod
    def tmp_files(cls, target_dir: Path, *names: str, mkdir=False):
        files = [target_dir / n for n in names]
        return cls(
            pass_ext=None,
            fail_ext=None,
            target_dir=target_dir,
            mkdir=mkdir,
            files=files,
        )

    def __enter__(self):
        if self._mkdir:
            self._target_dir.mkdir(exist_ok=False)

        for f in self._files:
            f.touch(exist_ok=False)

        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        for f in self._files:
            f.unlink(missing_ok=False)

        if self._mkdir:
            self._target_dir.rmdir()


def _check_options_generator(method, directory: Path, extension: str):
    with TemporaryFiles.tmp_in_dir(directory, extension) as c:
        options = method()
        assert c._pass_ext.name in options
        assert options[c._pass_ext.name] == c._pass_ext
        assert c._fail_ext.name not in options
        with pytest.raises(KeyError):
            assert options[c._fail_ext.name] == c._fail_ext


def test__select_from_dict():
    k, v = "foo", "bar"
    assert _select_from_dict(k, {k: v}) == v


def test_pipeline_options_genome_reference(mock_genome_dir):
    assembly, release, mock_dir = mock_genome_dir

    with TemporaryFiles.tmp_files(
        mock_dir,
        "something.dna.toplevel.fa",
        "something.dna.toplevel.chrom",
        mkdir=True,
    ):
        options = PipelineUIOptions.genome_reference()
        assert "{} - Ensembl release {}".format(assembly, release) in options


def test_pipeline_options_annotated_mutations():
    _check_options_generator(
        PipelineUIOptions.annotated_mutations,
        Directories.app_annotated_inputs(),
        "anno",
    )


def test_pipeline_options_immunopeptidome():
    _check_options_generator(
        PipelineUIOptions.immunopeptidome,
        Directories.app_immunopeptidomes(),
        "bed",
    )


def test_pipeline_options_substitution_method():
    for v in (192, 7):
        assert v in PipelineUIOptions.substitution_method().values()


def test_pipeline_options_coordinates():
    _check_options_generator(
        PipelineUIOptions.coordinates,
        Directories.app_coordinate_files(),
        "bed",
    )


def test_pipeline_processing_genome_reference(mock_genome_dir):
    assembly, release, mock_dir = mock_genome_dir

    selection_string = "{} - Ensembl release {}".format(assembly, release)

    output = PipelineUIProcessing.genome_reference(selection_string)

    assert isinstance(output, GenomePaths)
    assert output.fasta.name.endswith(".fa")
    assert output.sizes.name.endswith(".chrom")


def test_pipeline_processing_annotated_mutations():
    with TemporaryFiles.tmp_in_dir(
        Directories.app_annotated_inputs(), "anno"
    ) as c:
        assert (
            PipelineUIProcessing.annotated_mutations(c._pass_ext.name)
            == c._pass_ext
        )

        with pytest.raises(KeyError):
            PipelineUIProcessing.annotated_mutations(c._fail_ext.name)


def test_pipeline_processing_immunopeptidome():
    with TemporaryFiles.tmp_in_dir(
        Directories.app_immunopeptidomes(), "bed"
    ) as c:
        assert (
            PipelineUIProcessing.immunopeptidome(c._pass_ext.name)
            == c._pass_ext
        )

        with pytest.raises(KeyError):
            PipelineUIProcessing.immunopeptidome(c._fail_ext.name)


def test_pipeline_processing_substitution_method():
    assert PipelineUIProcessing.substitution_method("SSB192") == 192
    assert PipelineUIProcessing.substitution_method("SSB7") == 7

    with pytest.raises(KeyError):
        PipelineUIProcessing.substitution_method("SSB000")


def test_pipeline_processing_job_name():
    cache_key = "SOPRANO_CACHE"
    if cache_key in os.environ:
        _soprano_cache = os.environ[cache_key]
        del os.environ[cache_key]
    else:
        _soprano_cache = _SOPRANO_DEFAULT_CACHE.as_posix()
    try:
        mock_dir = Path("/just/for/pytest")
        os.environ[cache_key] = mock_dir.as_posix()
        job_name = "soprano"
        assert PipelineUIProcessing.job_name(job_name) == mock_dir / job_name
    finally:
        os.environ[cache_key] = _soprano_cache
