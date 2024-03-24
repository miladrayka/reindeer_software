"""Module provides functions for Proteo-Chemometrics Interaction
Fingerprints (PrtCmm IFPs) feature vector."""

from src.feature_generators.oic_dwic import FeatureGenerator
from .IFP import pro_lig_ifp


class PrtCmmIFP(FeatureGenerator):
    """PrtCmmIFP is a class to generate features
    for Proteo-Chemometrics Interaction Fingerprints (PrtCmm IFPs).
    """

    def __init__(
        self,
        pathfiles: str,
        filename: str,
        ligand_format: str,
        parameters: dict,
    ) -> None:
        """_summary_

        Parameters
        ----------
        pathfiles : str
            Path to protein-ligand complexes folders.
        filename : str
            Filename to save generated features.
        ligand_format : str
            Ligand format file. For now, only mol2 is acceptable.
        parameters : dict
            parameters of PrtCmm IFPs
        """
        super().__init__(pathfiles, filename, ligand_format)
        self.parameters = parameters

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
        try:
            pl_inter = pro_lig_ifp(
                fn_pro=protein_file,
                fn_lig=ligand_file,
                pdb="PDBID",
                int_cutoff=self.parameters["contact_para"]["int_cutoff"],
                ch_para=self.parameters["ch_para"],
                contact_bins=self.parameters["contact_para"]["bins"],
            )

            pl_inter.find_contacts()

            fv = pl_inter.featurize_ifp_mulpr(
                contact_bins=None,
                ifptype=self.parameters["ifp_para"]["ifptype"],
                alg_type=self.parameters["alg_type"],
                ecfp_radius=self.parameters["ifp_para"]["degrees"],
                heavy_atoms=self.parameters["ifp_para"]["heavy_atoms"],
                base_prop=self.parameters["ifp_para"]["prop"],
                hash_type=self.parameters["ifp_para"]["hash_type"],
                idf_power=self.parameters["ifp_para"]["idf_power"],
                folding_para=self.parameters["folding_para"],
            )

            feature_vector = fv[self.parameters["folding_para"]["power"][0]]
            feature_names = ["bit_" + str(item) for item in range(0, 1024)]

            return dict(zip(feature_names, feature_vector))

        except RecursionError:

            feature_vector = [0] * 1024
            feature_names = ["bit_" + str(item) for item in range(0, 1024)]
            return dict(zip(feature_names, feature_vector))
