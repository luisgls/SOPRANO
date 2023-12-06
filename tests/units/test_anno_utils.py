import pytest

from SOPRANO.utils import anno_utils


class TestFindVCFs:
    @staticmethod
    def test_vcf_gz_detected(tmp_path):
        vcf_path = tmp_path / "test.vcf.gz"
        vcf_path.touch()
        assert anno_utils.find_vcf_files(vcf_path) == [vcf_path]

    @staticmethod
    def test_vcf_detected(tmp_path):
        vcf_path = tmp_path / "test.vcf"
        vcf_path.touch()
        assert anno_utils.find_vcf_files(vcf_path) == [vcf_path]

    @staticmethod
    def test_not_vcf_gz_undetected(tmp_path):
        vcf_path = tmp_path / "test.foo.gz"
        vcf_path.touch()

        with pytest.raises(anno_utils.NotVCF):
            anno_utils.find_vcf_files(vcf_path)

    @staticmethod
    def test_multi_vcf_gz_detected(tmp_path):
        vcf_paths = [tmp_path / f"{tag}.vcf.gz" for tag in ("a", "b")]
        for p in vcf_paths:
            p.touch()

        expected = sorted(vcf_paths)
        found = sorted(anno_utils.find_vcf_files(tmp_path))

        assert expected == found

    @staticmethod
    def test_mixture_of_exts(tmp_path):
        unacceptable_paths = [
            tmp_path / "bad.tool",
            tmp_path / "partial.vcf",
            tmp_path / "woof.gz",
        ]
        acceptable_paths = [tmp_path / "ok.vcf.gz", tmp_path / "ok.VCF"]
        all_vcf_paths = unacceptable_paths + acceptable_paths

        for vcf_path in all_vcf_paths:
            vcf_path.touch()

        # order is not important, but may be shuffled by sets
        detected_paths = anno_utils.find_vcf_files(tmp_path).sort()
        expected_paths = acceptable_paths.sort()

        assert detected_paths == expected_paths

    @staticmethod
    def test_none_detected(tmp_path):
        with pytest.raises(anno_utils.NoVCFs):
            anno_utils.find_vcf_files(tmp_path)
