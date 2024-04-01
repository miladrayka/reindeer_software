"""Module provides functions for generating Multi Shell Occurence of
Interatomic Contacts feature vector from OnionNet-2 scoring function."""

from itertools import product
from collections import OrderedDict

import numpy as np
from scipy.spatial.distance import cdist
from Bio.PDB import PDBParser

from src.script import mol2parser
from .oic_dwic import FeatureGenerator


class MultiShellOIC(FeatureGenerator):
    """MultiShellOTC is a class to generate features for OnionNet-2 scoring
    function. This code is developed by using original github code
    (https://github.com/zchwang/OnionNet-2).
    """

    def __init__(
        self, pathfiles: str, filename: str, ligand_format: str, n_shells: int
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
        n_shells : int
            Number of shells.
        """
        super().__init__(pathfiles, filename, ligand_format)
        self.n_shells = n_shells
        self.defined_residues = [
            "GLY",
            "ALA",
            "VAL",
            "LEU",
            "ILE",
            "PRO",
            "PHE",
            "TYR",
            "TRP",
            "SER",
            "THR",
            "CYS",
            "MET",
            "ASN",
            "GLN",
            "ASP",
            "GLU",
            "LYS",
            "ARG",
            "HIS",
            "OTH",
        ]
        self.defined_elements = ["H", "C", "O", "N", "P", "S", "Hal", "DU"]
        self.keys = [
            "_".join(x)
            for x in list(product(self.defined_residues, self.defined_elements))
        ]
        self.feature_names = [
            key + "_" + str(shell)
            for shell in range(1, n_shells + 1)
            for key in self.keys
        ]

    def features_generator(self, ligand_file: str, protein_file: str) -> dict:
        """Generate feature for a protein-ligand complex.

        Parameters
        ----------
        ligand_file : str
            File of ligand in mol2.
        protein_file : str
            File of protein in pdb.

        Returns
        -------
        dict
            Return feature vector of a protein-ligand complex.
        """

        ligand_element_list, ligand_coords_list = self._loadmol2(ligand_file)
        residue_list, all_residue_coords_list = self._loadpdb(protein_file)
        residue_atom_dist, residue_atom_pairs = self._calculate_distance(
            residue_list,
            all_residue_coords_list,
            ligand_element_list,
            ligand_coords_list,
        )
        feature_vector = self._count_contacts(
            residue_atom_dist, residue_atom_pairs
        )

        return dict(zip(self.feature_names, feature_vector))

    def _loadmol2(self, ligand_file: str) -> tuple:
        """Load a ligand file in mol2.

        Parameters
        ----------
        ligand_file : str
            File of ligand in mol2.

        Returns
        -------
        tuple
            A tuple of list of all ligand atoms and their coordinates.
        """

        ligand = mol2parser.Mol2Parser(ligand_file)
        ligand.parse()

        ligand_element_list = list(
            map(lambda x: x[0], ligand.molecule_info["atom_name"].values())
        )

        for item, _ in enumerate(ligand_element_list):

            if ligand_element_list[item] in ["F", "Cl", "Br", "I"]:
                ligand_element_list[item] = "Hal"

            elif ligand_element_list[item] not in ["H", "C", "O", "N", "P", "S"]:
                ligand_element_list[item] = "DU"

            else:
                continue

        ligand_coords_list = np.array(
            list(ligand.molecule_info["coords"].values())
        ).astype(np.float32)

        return (ligand_element_list, ligand_coords_list * 0.1)

    def _loadpdb(self, protein_file: str) -> tuple:
        """Load a protein file in pdb.

        Parameters
        ----------
        protein_file : str
            File of protein in pdb.

        Returns
        -------
        tuple
            A tuple of list of all residues, their heavy atoms, and their
            coordinates.
        """

        parser = PDBParser(PERMISSIVE=True, QUIET=True)
        protein = parser.get_structure("", protein_file)

        residue_list = []
        all_residue_coords_list = []

        for residue in protein.get_residues():

            residue_coords = []
            if residue.get_resname() in self.defined_residues:
                residue_list.append(residue.get_resname())
            else:
                residue_list.append("OTH")

            for atom in residue.get_atoms():

                if atom.element != "H":
                    residue_coords.append(list(atom.get_coord() * 0.1))

            all_residue_coords_list.append(np.array(residue_coords))

        return (residue_list, all_residue_coords_list)

    def _calculate_distance(
        self,
        residue_list: list,
        all_residue_coords_list: np.array,
        ligand_element_list: list,
        ligand_coords_list: np.array,
    ) -> tuple:
        """Calculate distance between ligand and protein atoms.

        Parameters
        ----------
        residue_list : list
            A list of all protein residues.
        all_residue_coords_list : np.array
            A numpy array of all protein atoms and their corrdinates.
        ligand_element_list : list
            A list of all ligand atoms.
        ligand_coords_list : np.array
            A numpy array of all ligand coordinates.

        Returns
        -------
        tuple
            A tuple including a numpy array of minimum distances between
            residues and ligand atoms, and a list of strings of each
            residue-atom pair.
        """

        residue_atom_dist = []
        residue_atom_pairs = []

        for res, res_coords in zip(residue_list, all_residue_coords_list):
            for ele, atom_coords in zip(ligand_element_list, ligand_coords_list):
                pair = f"{res}_{ele}"
                dist_mtx = cdist(
                    atom_coords.reshape(1, -1), res_coords, metric="euclidean"
                )
                residue_atom_pairs.append(pair)
                residue_atom_dist.append(dist_mtx.min())

        residue_atom_dist = np.array(residue_atom_dist)

        return (residue_atom_dist, residue_atom_pairs)

    def _count_contacts(
        self,
        residue_atom_dist: np.array,
        residue_atom_pairs: list,
    ) -> np.array:
        """_summary_

        Parameters
        ----------
        residue_atom_dist : np.array
            A numpy array of minimum distances between residues and ligand
            atoms.
        residue_atom_pairs : list
            A list of strings of each residue-atom pair.

        Returns
        -------
        np.array
            Counted countact for each residue-atom pair in a specific shell.
        """

        outermost = 0.05 * (self.n_shells + 1)
        ncutoffs = np.linspace(0.1, outermost, self.n_shells)

        temp_counts = []
        onion_counts = []

        for i, cutoff in enumerate(ncutoffs):
            contact_bool = (residue_atom_dist <= cutoff) * 1
            if i == 0:
                onion_counts.append(contact_bool)
            else:
                onion_counts.append(contact_bool - temp_counts[-1])
            temp_counts.append(contact_bool)
        temp_counts = []

        results = []

        for n in range(len(ncutoffs)):
            d = OrderedDict()
            d = d.fromkeys(self.keys, 0)
            for e_e, c in zip(residue_atom_pairs, onion_counts[n]):
                d[e_e] += c
            results.append(np.array(list(d.values())).ravel())
        results = np.concatenate(results, axis=0)

        return results
