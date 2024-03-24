"""Module provides functions for generating Extended Connectivity Interaction Features
(ECIF) feature vector."""

from itertools import product

import pandas as pd
from rdkit import Chem
from scipy.spatial.distance import cdist
from Bio.PDB import PDBParser

from src.script import utils
from src.script import Mol
from .oic_dwic import FeatureGenerator


class ECIF(FeatureGenerator):
    """ECIF is a class to generate features for ECIF techniques.
    This code is developed by using ECIF original github code (https://github.com/DIFACQUIM/ECIF)
    """

    def __init__(
        self, pathfiles: str, filename: str, ligand_format: str, cutoff: float
    ) -> None:
        """
        Parameters
        ----------
        pathfiles : str
            Path to protein-ligand complexes folders.
        filename : str
            Filename to save generated features.
        ligand_format : str
            Ligand format file. For now, only sdf is acceptable.
        cutoff : float
            Distance threshold.
        """
        super().__init__(pathfiles, filename, ligand_format)
        self.cutoff = cutoff
        self.possible_ecif = [
            i[0] + "-" + i[1]
            for i in product(utils.ECIF_ProteinAtoms, utils.ECIF_LigandAtoms)
        ]
        self.atom_keys = pd.read_csv("./src/files/PDB_Atom_Keys.csv", sep=",")

    def features_generator(self, ligand_file: str, protein_file: str) -> dict:
        """Generate feature for a protein-ligand complex.

        Parameters
        ----------
        ligand_file : str
            File of ligand in sdf.
        protein_file : str
            File of protein in pdb.

        Returns
        -------
        dict
            Return feature vector of a protein-ligand complex.
        """

        pairs = self._getplpairs(protein_file, ligand_file, self.cutoff)
        feature_vector = [list(pairs["ECIF_PAIR"]).count(x) for x in self.possible_ecif]
        return dict(zip(self.possible_ecif, feature_vector))

    def _getplpairs(self, protein_file: str, ligand_file: str, cutoff: float) -> list:
        """Count pairs of possible ecif atom types.

        Parameters
        ----------
        protein_file : str
            File of protein in pdb.
        ligand_file : str
            File of ligand in sdf.
        cutoff : float
            Distance threshold.

        Returns
        -------
        list
            Counted number of possible ecif pairs.
        """

        target = self._loadpdb(protein_file)
        ligand = self._loadsdf(ligand_file)

        pairs = list(product(target["ECIF_ATOM_TYPE"], ligand["ECIF_ATOM_TYPE"]))
        pairs = [x[0] + "-" + x[1] for x in pairs]
        pairs = pd.DataFrame(pairs, columns=["ECIF_PAIR"])
        distances = cdist(
            target[["X", "Y", "Z"]], ligand[["X", "Y", "Z"]], metric="euclidean"
        ).reshape(-1, 1)
        pairs["DISTANCE"] = distances

        pairs = pairs[pairs["DISTANCE"] <= cutoff].reset_index(drop=True)

        return pairs

    def _loadpdb(self, protein_file: str) -> pd.DataFrame:
        """Load PDB file.
        Parameters
        ----------
        protein_file : str
            File of protein in pdb.

        Returns
        -------
        pd.DataFrame
            Return some information of a protein in a dataframe.
        """

        parser = PDBParser(PERMISSIVE=True, QUIET=True)
        protein = parser.get_structure("", protein_file)
        ecif_atoms = []
        atom_index = 1

        for residue in protein.get_residues():
            if residue.id[0] == " ":
                for atom in residue.get_atoms():
                    if atom.element in ["C", "N", "O", "S"]:
                        ecif_atoms.append(
                            [
                                atom_index,
                                residue.get_resname() + "-" + atom.name,
                                atom.coord[0],
                                atom.coord[1],
                                atom.coord[2],
                            ]
                        )
                        atom_index += 1

        df = pd.DataFrame(ecif_atoms, columns=["ATOM_INDEX", "PDB_ATOM", "X", "Y", "Z"])
        df = (
            df.merge(self.atom_keys, left_on="PDB_ATOM", right_on="PDB_ATOM")[
                ["ATOM_INDEX", "ECIF_ATOM_TYPE", "X", "Y", "Z"]
            ]
            .sort_values(by="ATOM_INDEX")
            .reset_index(drop=True)
        )
        df.iloc[:, 2:] = df.iloc[:, 2:].astype(float).round(3)

        return df

    def _loadsdf(self, ligand_file: str) -> pd.DataFrame:
        """Load SDF file.

        Parameters
        ----------
        ligand_file : str
            File of ligand in sdf.

        Returns
        -------
        pd.DataFrame
            Return some information of a ligand in a dataframe.
        """

        m = Chem.MolFromMolFile(ligand_file, sanitize=False)
        m.UpdatePropertyCache(strict=False)

        ecif_atoms = []

        for atom in m.GetAtoms():
            if atom.GetSymbol() != "H" and atom.GetSymbol() in [
                "C",
                "N",
                "O",
                "S",
                "P",
                "F",
                "Cl",
                "Br",
                "I",
            ]:
                entry = [int(atom.GetIdx())]
                entry.append(self._getatom_type(atom))
                pos = m.GetConformer().GetAtomPosition(atom.GetIdx())
                entry.append(float(f"{pos.x:.4f}"))
                entry.append(float(f"{pos.y:.4f}"))
                entry.append(float(f"{pos.z:.4f}"))
                ecif_atoms.append(entry)

        df = pd.DataFrame(ecif_atoms)
        df.columns = ["ATOM_INDEX", "ECIF_ATOM_TYPE", "X", "Y", "Z"]

        return df

    def _getatom_type(self, atom: Mol) -> str:
        """Generate properties for a given atom in RDKit Mol object.

        Parameters
        ----------
        atom : Mol
            Atom RDKit Mol object.

        Returns
        -------
        str
            A string of properties of an atom.
        """

        atom_type = [
            atom.GetSymbol(),
            str(atom.GetExplicitValence()),
            str(
                len(
                    [x.GetSymbol() for x in atom.GetNeighbors() if x.GetSymbol() != "H"]
                )
            ),
            str(
                len(
                    [x.GetSymbol() for x in atom.GetNeighbors() if x.GetSymbol() == "H"]
                )
            ),
            str(int(atom.GetIsAromatic())),
            str(int(atom.IsInRing())),
        ]

        return ";".join(atom_type)
