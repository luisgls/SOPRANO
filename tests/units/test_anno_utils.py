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
