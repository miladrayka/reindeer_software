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

# Possible protein atom types for ECIF.
ECIF_ProteinAtoms = [
    "C;4;1;3;0;0",
    "C;4;2;1;1;1",
    "C;4;2;2;0;0",
    "C;4;2;2;0;1",
    "C;4;3;0;0;0",
    "C;4;3;0;1;1",
    "C;4;3;1;0;0",
    "C;4;3;1;0;1",
    "C;5;3;0;0;0",
    "C;6;3;0;0;0",
    "N;3;1;2;0;0",
    "N;3;2;0;1;1",
    "N;3;2;1;0;0",
    "N;3;2;1;1;1",
    "N;3;3;0;0;1",
    "N;4;1;2;0;0",
    "N;4;1;3;0;0",
    "N;4;2;1;0;0",
    "O;2;1;0;0;0",
    "O;2;1;1;0;0",
    "S;2;1;1;0;0",
    "S;2;2;0;0;0",
]

# Possible ligand atom types for ECIF.
ECIF_LigandAtoms = [
    "Br;1;1;0;0;0",
    "C;3;3;0;1;1",
    "C;4;1;1;0;0",
    "C;4;1;2;0;0",
    "C;4;1;3;0;0",
    "C;4;2;0;0;0",
    "C;4;2;1;0;0",
    "C;4;2;1;0;1",
    "C;4;2;1;1;1",
    "C;4;2;2;0;0",
    "C;4;2;2;0;1",
    "C;4;3;0;0;0",
    "C;4;3;0;0;1",
    "C;4;3;0;1;1",
    "C;4;3;1;0;0",
    "C;4;3;1;0;1",
    "C;4;4;0;0;0",
    "C;4;4;0;0;1",
    "C;5;3;0;0;0",
    "C;5;3;0;1;1",
    "C;6;3;0;0;0",
    "Cl;1;1;0;0;0",
    "F;1;1;0;0;0",
    "I;1;1;0;0;0",
    "N;3;1;0;0;0",
    "N;3;1;1;0;0",
    "N;3;1;2;0;0",
    "N;3;2;0;0;0",
    "N;3;2;0;0;1",
    "N;3;2;0;1;1",
    "N;3;2;1;0;0",
    "N;3;2;1;0;1",
    "N;3;2;1;1;1",
    "N;3;3;0;0;0",
    "N;3;3;0;0;1",
    "N;3;3;0;1;1",
    "N;4;1;2;0;0",
    "N;4;1;3;0;0",
    "N;4;2;1;0;0",
    "N;4;2;2;0;0",
    "N;4;2;2;0;1",
    "N;4;3;0;0;0",
    "N;4;3;0;0;1",
    "N;4;3;1;0;0",
    "N;4;3;1;0;1",
    "N;4;4;0;0;0",
    "N;4;4;0;0;1",
    "N;5;2;0;0;0",
    "N;5;3;0;0;0",
    "N;5;3;0;1;1",
    "O;2;1;0;0;0",
    "O;2;1;1;0;0",
    "O;2;2;0;0;0",
    "O;2;2;0;0;1",
    "O;2;2;0;1;1",
    "P;5;4;0;0;0",
    "P;6;4;0;0;0",
    "P;6;4;0;0;1",
    "P;7;4;0;0;0",
    "S;2;1;0;0;0",
    "S;2;1;1;0;0",
    "S;2;2;0;0;0",
    "S;2;2;0;0;1",
    "S;2;2;0;1;1",
    "S;3;3;0;0;0",
    "S;3;3;0;0;1",
    "S;4;3;0;0;0",
    "S;6;4;0;0;0",
    "S;6;4;0;0;1",
    "S;7;4;0;0;0",
]


# Parameters for PrtCmm IFP and PLEC FP:
# alg_type - ifp algorithm ('avg', 'classic')
# int_cutoff: distance cutoff for identifying contacting atoms in protein-ligand interfaces
# ch_para - parameters for constructing alpha shapes (concave hulls) for protein and ligand
#                weighted: weighted alpha shape (1) or mlnot (0)
#                alpha: alpha values for protein (alpha[0]) and ligand (alpha[1])
#                alpha_step: steps for tuning alpha shapes if alpha = -1 (alpha_step[0] for protein and alpha_step[1] fro ligand)
# bins - for finding protein-ligand contacts in different ranges
# ifptype - interaction fingerprint type, 'splif', 'ecfp' or 'plec'
# degrees - ECFP radii for the ifp
# prop - a list of atomic properties, full list as below
# heavy_atoms: use heavy atoms or all atoms
# hash_type: type for the hash function ('str' or 'vec')
# idf_power: power for the identifiers hashed by hash_ecfp ('str')
# folder_para - parameters for fingerprint folding
#                power: fingerprint folding size (2^power bits)
#                counts: use occurences of identifiers (1) or not (0)

prtcmm_ifp_parameters = {
    "alg_type": "avg",
    "ch_para": {"weighted": 0, "alpha": [-1, -1], "alpha_step": [0.1, 0.1]},
    "contact_para": {"int_cutoff": 4.5, "bins": [(0, 4.5)]},
    "ifp_para": {
        "ifptype": "ecfp",
        "degrees": [1, 1],
        "prop": [
            "AtomicMass",
            "TotalConnections",
            "HeavyNeighborCount",
            "HCount",
            "FormalCharge",
        ],
        "heavy_atoms": 1,
        "hash_type": "vec",
        "idf_power": 64,
    },
    "folding_para": {"power": [9], "counts": 1},
}

plec_fp_parameters = {
    "alg_type": "avg",
    "ch_para": {"weighted": 0, "alpha": [-1, -1], "alpha_step": [0.1, 0.1]},
    "contact_para": {"int_cutoff": 4.5, "bins": [(0, 4.5)]},
    "ifp_para": {
        "ifptype": "plec",
        "degrees": [5, 1],
        "prop": [
            "AtomicMass",
            "TotalConnections",
            "HeavyNeighborCount",
            "HCount",
            "FormalCharge",
        ],
        "heavy_atoms": 1,
        "hash_type": "vec",
        "idf_power": 64,
    },
    "folding_para": {"power": [11], "counts": 1},
}

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
