"""Command Line Interface (CLI) for REINDEERS software."""

import argparse

from reindeer.feature_generators import oic_dwic
from reindeer.feature_generators import ecif
from reindeer.feature_generators import ms_oic
from reindeer.script import utils

if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="""Generate features for
                                     set of given structures"""
    )
    parser.add_argument(
        "-m",
        "--method",
        type=str,
        help="""Feature generation method.
                        Only OIC, DWIC, ECIF, and MS-OIC
                        are implemented for now.""",
        required=True,
    )
    parser.add_argument(
        "-d",
        "--directory",
        type=str,
        help="""directory of
                        structures files""",
        required=True,
    )
    parser.add_argument(
        "-f",
        "--file_name",
        type=str,
        help="""Name for saving
                        generated features.""",
        required=True,
    )
    parser.add_argument(
        "-n",
        "--n_jobs",
        type=int,
        help="""Number of cpu cores
                        for parallelization""",
        required=True,
    )

    args = parser.parse_args()

    if args.method == "OIC":

        oic = oic_dwic.InterAtomicContact(
            pathfiles=args.directory,
            filename=args.file_name,
            ligand_format="mol2",
            amino_acid_classes=utils.amino_acid_classes_OIC,
            cutoff=12.0,
            feature_type="OIC",
            exp=None,
        )
        oic.generate_features(n_jobs=args.n_jobs)

    elif args.method == "DWIC":

        dwic = oic_dwic.InterAtomicContact(
            pathfiles=args.directory,
            filename=args.file_name,
            ligand_format="mol2",
            amino_acid_classes=utils.amino_acid_classes_DWIC,
            cutoff=12.0,
            feature_type="DWIC",
            exp=2,
        )
        dwic.generate_features(n_jobs=args.n_jobs)

    elif args.method == "ECIF":
        ecif = ecif.ECIF(
            pathfiles=args.directory,
            filename=args.file_name,
            ligand_format="sdf",
            cutoff=6.0,
        )
        ecif.generate_features(n_jobs=args.n_jobs)

    elif args.method == "MS-OIC":
        ms_oic = ms_oic.MultiShellOIC(
            pathfiles=args.directory,
            filename=args.file_name,
            ligand_format="mol2",
            n_shells=62,
        )
        ms_oic.generate_features(n_jobs=args.n_jobs)

    else:
        raise ValueError("Method should be selected from implemented methods.")
