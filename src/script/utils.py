"""Module provides some useful functions."""

from collections import defaultdict
from Bio.PDB import PDBParser
from tqdm.auto import tqdm
from joblib import Parallel

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


class ProgressParallel(Parallel):
    """A wrapper for Parallel class of joblib. This wrapper
    add progress bar to Parallel. Code source: https://stackoverflow.com/
    a/61900501.
    """
    def __init__(self, *args, use_tqdm=True, total=None, **kwargs):
        self._use_tqdm = use_tqdm
        self._total = total
        self._pbar = None
        super().__init__(*args, **kwargs)

    def __call__(self, *args, **kwargs):
        with tqdm(disable=not self._use_tqdm, total=self._total) as self._pbar:
            return Parallel.__call__(self, *args, **kwargs)

    def print_progress(self):
        if self._total is None:
            self._pbar.total = self.n_dispatched_tasks
        self._pbar.n = self.n_completed_tasks
        self._pbar.refresh()
