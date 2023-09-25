import pathlib

import streamlit as st

# Find genome options
_HOMO_SAPIENS_DIR = pathlib.Path(__file__).parent / "data" / "homo_sapiens"

_GENOME_DIRS = [item for item in _HOMO_SAPIENS_DIR.glob("*") if item.is_dir()]

_GENOME_NAMES = [
    "{} - Ensembl release {}".format(*x.name.split("_")[::-1])
    for x in _GENOME_DIRS
]

_GENOME_DICT = {
    name: dir_path for name, dir_path in zip(_GENOME_NAMES, _GENOME_DIRS)
}

_ANNO_DIR = pathlib.Path(__file__).parent / "examples"
_ANNO_OPTIONS = {x.name: x for x in _ANNO_DIR.glob("*.anno*")}

_BED_DIR = pathlib.Path(__file__).parent / "immunopeptidomes" / "human"
_BED_OPTIONS = {x.name: x for x in _BED_DIR.glob("*.bed")}

# Init app
st.title("SOPRANO")

# Derived genome definitions
genome_selection = st.selectbox(
    "Select a reference genome:", _GENOME_DICT.keys()
)
ref_id, rel_id = genome_selection.split(" - Ensembl release ")
toplevel_path = (
    _GENOME_DICT[genome_selection] / f"Homo_sapiens.{ref_id}.dna.toplevel.fa"
)
st.text(f"Selected: {toplevel_path}")

# VEP annotated file
input_selection = st.selectbox(
    "Select a VEP annotated file:", _ANNO_OPTIONS.keys()
)
input_path = _ANNO_OPTIONS[input_selection]
st.text(f"Selected: {input_path}")

# BED protein transcript file
bed_selection = st.selectbox("Select a BED protein file:", _BED_OPTIONS.keys())
bed_path = _BED_OPTIONS[bed_selection]
st.text(f"Selected: {bed_path}")
