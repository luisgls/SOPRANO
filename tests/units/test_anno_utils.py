import pytest

from SOPRANO.utils import anno_utils


class TestFindVCFs:
    @staticmethod
    def test_single_file_ok(tmp_path):
        vcf_path = tmp_path / "test.vcf.gz"
        vcf_path.touch()
        assert anno_utils.find_vcf_files(vcf_path) == [vcf_path]

    @staticmethod
    def test_single_file_no_gz(tmp_path):
        vcf_path = tmp_path / "test.vcf.bar"
        vcf_path.touch()

        with pytest.raises(anno_utils.NotGZ):
            anno_utils.find_vcf_files(vcf_path)

    @staticmethod
    def test_single_file_no_vcf(tmp_path):
        vcf_path = tmp_path / "test.foo.gz"
        vcf_path.touch()

        with pytest.raises(anno_utils.NotVCF):
            anno_utils.find_vcf_files(vcf_path)

    @staticmethod
    def test_multi_file_ok(tmp_path):
        vcf_paths = [tmp_path / f"{tag}.vcf.gz" for tag in ("a", "b")]
        for p in vcf_paths:
            p.touch()

        expected = sorted(vcf_paths)
        found = sorted(anno_utils.find_vcf_files(tmp_path))

        assert expected == found

    @staticmethod
    def test_some_files_ok(tmp_path):
        bad_paths = [
            tmp_path / "bad.tool",
            tmp_path / "partial.vcf",
            tmp_path / "woof.gz",
        ]
        good_paths = [tmp_path / "ok.vcf.gz", tmp_path / "ok.VCF.Gz"]
        vcf_paths = bad_paths + good_paths

        for p in vcf_paths:
            p.touch()

        # order is not important, but contets is
        found = anno_utils.find_vcf_files(tmp_path).sort()
        expected = good_paths.sort()

        assert found == expected

    @staticmethod
    def test_multi_file_fail(tmp_path):
        with pytest.raises(anno_utils.NoVCFs):
            anno_utils.find_vcf_files(tmp_path)
