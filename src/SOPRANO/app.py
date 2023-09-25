import pathlib

import streamlit as st

# Find genome options
homo_sapiens_dir = pathlib.Path(__file__).parent / "data" / "homo_sapiens"

_GENOME_DIRS = [item for item in homo_sapiens_dir.glob("*") if item.is_dir()]

_GENOME_NAMES = [
    "{} - Ensembl release {}".format(*x.name.split("_")[::-1])
    for x in _GENOME_DIRS
]

_GENOME_DICT = {
    name: dir_path for name, dir_path in zip(_GENOME_NAMES, _GENOME_DIRS)
}

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
