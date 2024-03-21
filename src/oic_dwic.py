"""Module provides functions for generating DWIC and OIC feature vector."""

import time
from glob import glob
from pathlib import Path
from joblib import Parallel, delayed

import numpy as np
import pandas as pd
from scipy import spatial

from script import utils
from script import mol2parser


class FeatureGenerator(object):
    """Generate feature for protein-ligand complexes."""

    def __init__(self, pathfiles: str, filename: str, ligand_format: str) -> None:
        """_summary_

        Parameters
        ----------
        pathfiles : str
            Path to protein-ligand complexes folders.
        filename : str
            Filename to save generated features.
        ligand_format : str
            Ligand format file. For now, only mol2 is acceptable.
        """

        self.pathfiles = pathfiles
        self.ligand_format = ligand_format
        self.filename = filename

    def generate_features(self, n_jobs: int = 1) -> None:
        """Generate features parallelly and save them to filename.csv.

        Parameters
        ----------
        n_jobs : int, default: 1
            The maximum number of concurrently running jobs.  If -1 all CPUs are used.
        """

        entries = Path(self.pathfiles)
        start = time.time()

        fvs = Parallel(n_jobs)(
            delayed(self.features_generator)(
                glob(str(entry) + f"\\*.{self.ligand_format}")[0],
                glob(str(entry) + "\\*.pdb")[0],
            )
            for entry in entries.iterdir()
        )

        names = [entry.name for entry in entries.iterdir()]

        end = time.time()
        m, s = divmod(end - start, 60)
        h, m = divmod(m, 60)
        print(f"Time: {h}:{m}:{s:.2f}")

        self._save_to_csv(dict(zip(names, fvs)), self.filename)

    def _save_to_csv(self, feature_vectors: dict, filename: str) -> None:
        """Save generated features to filename.csv

        Parameters
        ----------
        feature_vectors : dict
            Generated features.
        filename : str
            Filename to save generated features.
        """

        df = pd.DataFrame.from_dict(feature_vectors, orient="index")
        df = df.fillna(0.0)

        df.to_csv(filename, index=True, index_label="pdbid")

    def features_generator(self, ligand_file: str, protein_file: str) -> dict:
        """Feature generator function. Override this for different feature vector generator.

        Parameters
        ----------
        ligand_file : str
            File of ligand in mol2.
        protein_file : str
            File of protein in pdb.
        """


class InterAtomicContact(FeatureGenerator):
    """InterAtomicContact is a class to generate features
    for both occurrence of interatomic contacts (OIC) and
    distance-weighted interatomic contacts (DWIC) techniques.
    """

    def __init__(
        self,
        pathfiles: str,
        filename: str,
        ligand_format: str,
        amino_acid_classes: list,
        cutoff: int,
        feature_type: str,
        exp: int,
    ) -> None:
        """
        Parameters
        ----------
        pathfiles : str
            Path to protein-ligand complexes folders.
        filename : str
            Filename to save generated features.
        ligand_format : str
            Ligand format file. For now, only mol2 is acceptable.
        amino_acid_classes : list
            List (or groups) of amino acids.
        cutoff : int
            Distance threshold.
        feature_type : str
            Generate features based on OIC or DWIC.
        exp : int
            Exponent variable in the DWIC equation.
        """

        super().__init__(pathfiles, filename, ligand_format)
        self.amino_acid_classes = amino_acid_classes
        self.cutoff = cutoff
        self.feature_type = feature_type
        self.exp = exp

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

        ligand = mol2parser.Mol2Parser(ligand_file)
        ligand.parse()

        extracted_protein_atoms = utils.pdb_parser(
            protein_file, self.amino_acid_classes
        )

        groups = extracted_protein_atoms.keys()
        feature_vector = {}

        for group in groups:
            for protein_element in utils.protein_elements:
                for ligand_element in utils.ligand_elements:

                    try:
                        ligand_coords = np.array(
                            ligand.get_molecule(ligand_element, "coords"), dtype="float"
                        ).reshape(-1, 3)

                    except TypeError:
                        continue

                    if extracted_protein_atoms[group][protein_element]:
                        protein_coords = np.array(
                            extracted_protein_atoms[group][protein_element],
                            dtype="float",
                        )

                        distances = spatial.distance.cdist(
                            ligand_coords, protein_coords
                        ).ravel()

                        if self.feature_type == "DWIC":
                            feature = self._weighting_sum(
                                distances, self.cutoff, self.exp
                            )

                        if self.feature_type == "OIC":
                            feature = self._counting_sum(distances, self.cutoff)

                        if len(groups) > 1:
                            name = group + "_" + protein_element + "_" + ligand_element

                        else:
                            name = protein_element + "_" + ligand_element

                        feature_vector[name] = feature

        return feature_vector

    def _weighting_sum(self, distances: np.array, cutoff: int, exp: int) -> float:
        """Method for generating DWIC feature.

        Parameters
        ----------
        distances : numpy array
            Calculated distances between protein and ligand atoms.
        cutoff : int
            Distance threshold.
        exp : int
            Exponent variable in the DWIC equation.

        Returns
        -------
        float
            Return a number as a feature.
        """

        selected_distances = distances[distances < cutoff]

        feature = sum(list(map(lambda x: 1.0 / (x**exp), selected_distances)))

        return feature

    def _counting_sum(self, distances: np.array, cutoff: int) -> float:
        """Method for generating OIC feature.

        Parameters
        ----------
        distances : numpy array
            Calculated distances between protein and ligand atoms.
        cutoff : int
            Distance threshold.

        Returns
        -------
        float
            Return a number as a feature.
        """

        selected_distances = distances[distances < cutoff]

        feature = len(selected_distances)

        return feature
