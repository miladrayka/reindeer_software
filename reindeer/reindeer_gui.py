"""This module provides a GUI for REINDEER software."""

import streamlit as st

from feature_generators import oic_dwic
from feature_generators import ecif
from feature_generators import ms_oic
from script import utils

st.set_page_config(
    page_title="REINDEER",
    page_icon=r".\logo\title_logo.png",
    layout="centered",
    initial_sidebar_state="collapsed",
)

st.title("REINDEER Software")
st.header("Protein-Ligand Feature Generator")
st.write(
    """**REINDEER** is a software for generating a feature vector
for a protein-ligand complex."""
)
st.image(r".\logo\Logo.png")

st.sidebar.header("Developer")
st.sidebar.write(
    """[GitHub](https://github.com/miladrayka/reindeer_software),
Developed by *[Milad Rayka](https://scholar.google.com/citations?user=NxF2f0cAAAAJ&hl=en)*."""
)
st.sidebar.divider()
st.sidebar.header("Citation")
st.sidebar.write(
    """**Reference**:
Paper is *under production.*"""
)

with st.expander("**Cautions**"):
    st.info("""1- Provided protein-ligand complex should have hydrogen atoms.""")
    st.info(
        """2- File formarts for protein and ligand are *.pdb* and *.mol2*.
    In the case of ECIF, instead of *.mol2*, *.sdf* file should be provided."""
    )
    st.info(
        """3- All protein-ligand complexes should be provided as the below example:\n
    ./test
    ├── 1a1e
    │   ├── 1a1e_ligand.mol2
    │   ├── 1a1e_ligand.sdf
    │   └── 1a1e_protein.pdb
    ├── 1a28
    │   ├── 1a28_ligand.mol2
    │   ├── 1a28_ligand.sdf
    │   └── 1a28_protein.pdb
    ├── 1a30
        ├── 1a30_ligand.mol2
        ├── 1a30_ligand.sdf
        └── 1a30_protein.pdb"""
    )

with st.expander("**Implemented Methods**"):
    st.write(
        """Current version of REINDEER provides four feature generation methods:"""
    )
    st.write(
        """1- Occurrence of Interatomic Contact (OIC) -
    **[Ref](https://academic.oup.com/bioinformatics/article/26/9/1169/199938?login=false)**"""
    )
    st.write(
        """2- Distance-Weighted Interatomic Contact (DWIC) -
    **[Ref](https://onlinelibrary.wiley.com/doi/abs/10.1002/minf.202060084)**"""
    )
    st.write(
        """3- Extended Connectivity Interaction Feature (ECIF) -
    **[Ref](https://academic.oup.com/bioinformatics/article/37/10/1376/5998664?login=false)**"""
    )
    st.write(
        """4- Multi-Shell Occurrence of Interatomic Contact (MS-OIC) -
    **[Ref](https://www.frontiersin.org/articles/10.3389/fchem.2021.753002/full)**"""
    )

path = st.text_input(
    "Enter the path of the folder of structures:",
    placeholder="Type a path...",
    help="Provide the path of protein-ligand structures, e.g., ./test/",
)
n_jobs = st.number_input(
    "Enter number of CPU core for parallelization:",
    min_value=-1,
    step=1,
    help="-1 value uses all CPU cores.",
)
filename = st.text_input(
    "Enter a name path for saving generated features: ",
    placeholder="Type a name...",
    help="Provide a name. Output file is in .csv format. For example, ./files/example",
)

st.write("Push the button for your desired feature vector:")

oic_condition = st.button("Generate OIC")

if oic_condition:
    with st.spinner("Please wait..."):
        oic = oic_dwic.InterAtomicContact(
            pathfiles=path,
            filename=f"{filename}_oic.csv",
            ligand_format="mol2",
            amino_acid_classes=utils.amino_acid_classes_OIC,
            cutoff=12.0,
            feature_type="OIC",
            exp=None,
        )
        oic.generate_features(n_jobs=n_jobs)
    st.success(f" Generated features are save into {filename}_oic.csv")

dwic_condition = st.button("Generate DWIC")
if dwic_condition:
    with st.spinner("Please wait..."):
        dwic = oic_dwic.InterAtomicContact(
            pathfiles=path,
            filename=f"{filename}_dwic.csv",
            ligand_format="mol2",
            amino_acid_classes=utils.amino_acid_classes_DWIC,
            cutoff=12.0,
            feature_type="DWIC",
            exp=2,
        )
        dwic.generate_features(n_jobs=n_jobs)
    st.success(f" Generated features are save into {filename}_dwic.csv")

ecif_condition = st.button("Generate ECIF")
if ecif_condition:
    with st.spinner("Please wait..."):
        ecif = ecif.ECIF(
            pathfiles=path,
            filename=f"{filename}_ecif.csv",
            ligand_format="sdf",
            cutoff=6.0,
        )
        ecif.generate_features(n_jobs=n_jobs)
    st.success(f" Generated features are save into {filename}_ecif.csv")

ms_oic_condition = st.button("Generate MS-OIC")
if ms_oic_condition:
    with st.spinner("Please wait..."):
        ms_oic = ms_oic.MultiShellOIC(
            pathfiles=path,
            filename=f"{filename}_ms_oic.csv",
            ligand_format="mol2",
            n_shells=62,
        )
        ms_oic.generate_features(n_jobs=n_jobs)
    st.success(f" Generated features are save into {filename}_ms_oic.csv")
