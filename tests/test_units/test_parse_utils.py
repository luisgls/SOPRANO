from argparse import Namespace
from unittest.mock import patch

from SOPRANO.utils.parse_utils import fix_species_arg, parse_genome_args


def test__fix_species_arg():
    ns = Namespace(species="Homo Sapiens")
    assert fix_species_arg(ns).species == "homo_sapiens"


def test_parse_genome_args(capsys):
    with patch("sys.argv", ["parse_genome_args"]):
        args = parse_genome_args()
        assert args.species == "homo_sapiens"
        assert args.assembly == "GRCh38"
        assert args.release == "110"
        assert args.primary_assembly is False
        assert args.download_only is False

    with patch(
        "sys.argv", ["parse_genome_args", "-s", "foo", "-a", "bar", "-p"]
    ):
        args = parse_genome_args()
        assert args.species == "foo"
        assert args.assembly == "bar"
        assert args.release == "110"
        assert args.primary_assembly is True
        assert args.download_only is False
