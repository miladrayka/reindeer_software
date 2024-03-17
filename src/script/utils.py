"""Module provides some useful functions."""

from collections import defaultdict
from Bio.PDB import PDBParser

# Amino acid list for Occurrence of interatomic contacts feature.
amino_acid_classes_OIC = [
    [
        "ARG",
        "LYS",
        "ASP",
        "GLU",
        "GLN",
        "ASN",
        "HIS",
        "SER",
        "THR",
        "CYS",
        "TRP",
        "TYR",
        "MET",
        "ILE",
        "LEU",
        "PHE",
        "VAL",
        "PRO",
        "GLY",
        "ALA",
    ]
]
# Amino acid list for Distance-weighted interatomic contacts features.
amino_acid_classes_DWIC = [
    ["ARG", "LYS", "ASP", "GLU"],
    ["GLN", "ASN", "HIS", "SER", "THR", "CYS"],
    ["TRP", "TYR", "MET"],
    ["ILE", "LEU", "PHE", "VAL", "PRO", "GLY", "ALA"],
]
# Ligand elements
ligand_elements = ["H", "C", "N", "O", "S", "P", "F", "Cl", "Br", "I"]

# Protein elements
protein_elements = ["H", "C", "N", "O", "S"]

# A pdb parser, writtein by BioPython, to extract protein atoms and
# coordinates.


def pdb_parser(pdbfile: str, amino_acid_classes: list) -> dict:
    """Extract protein atoms and their coordinates.

    Parameters
    ----------
        pdffile: str
            Protein file in .pdb format.
        amino_acid_classes: list
            List (or groups) of amino acids.

    Returns
    -------
        extracted_protein_atoms: dict
            Protein atoms and their coordinates in a dictionary.
    """
    number_of_groups = len(list(amino_acid_classes))

    amino_acid_groups = zip(
        [f"group{i}" for i in range(1, number_of_groups + 1)], amino_acid_classes
    )

    parser = PDBParser(PERMISSIVE=True, QUIET=True)
    protein = parser.get_structure("", pdbfile)

    extracted_protein_atoms = {
        f"group{i}": defaultdict(list) for i in range(1, number_of_groups + 1)
    }

    for name, group in amino_acid_groups:
        for residue in protein.get_residues():
            if residue.get_resname() in group:
                for atom in residue.get_atoms():
                    extracted_protein_atoms[name][atom.element].append(list(atom.coord))

    return extracted_protein_atoms
