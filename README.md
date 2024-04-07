[![Python 3.9](https://img.shields.io/badge/python-3.9-blue.svg)](https://www.python.org/downloads/release/python-360/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

# REINDEER Software
REINDEER is a software for the sturcture-based protein-ligand feature generation.

![Logo](https://github.com/miladrayka/reindeer_software/blob/main/reindeer/logo/Logo.png)

Currently, REINDEER provides only four below feature vector:

1- Occurrence of Interatomic Contact (OIC) - **[Ref](https://academic.oup.com/bioinformatics/article/26/9/1169/199938?login=false)**

2- Distance-Weighted Interatomic Contact (DWIC) - **[Ref](https://onlinelibrary.wiley.com/doi/abs/10.1002/minf.202060084)**

3- Extended Connectivity Interaction Feature (ECIF) - **[Ref](https://academic.oup.com/bioinformatics/article/37/10/1376/5998664?login=false)**

4- Multi-Shell Occurrence of Interatomic Contact (MS-OIC) - **[Ref](https://www.frontiersin.org/articles/10.3389/fchem.2021.753002/full)**

## Citation
Paper is under construction.

## Contact
Milad Rayka, milad.rayka@yahoo.com

## Notes

1- Provided protein-ligand complex should have hydrogen atoms

2- File formarts for protein and ligand are *.pdb* and *.mol2*. 
In the case of ECIF, instead of *.mol2*, *.sdf* file should be provided.

3- All protein-ligand complexes should be provided as the below example:

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
        └── 1a30_protein.pd


## Graphical User Interface (GUI)
REINDEER provides a GUI to make for feature generation methods. Check the [Tutorial](https://github.com/miladrayka/reindeer_software/blob/main/Tutorial.pdf) file.

![GUI](https://github.com/miladrayka/reindeer_software/blob/main/GUI_img.PNG))

## Command Line Interface (CLI)
REINDEER provides a simple CLI to ease the usage.

```
usage: reindeer_software.py [-h] -m METHOD -d DIRECTORY -f FILE_NAME -n N_JOBS

Generate features for set of given structures

optional arguments:
  -h, --help            show this help message and exit
  -m METHOD, --method METHOD
                        Feature generation method. Only OIC, DWIC, ECIF, and
                        MS-OIC are implemented for now.
  -d DIRECTORY, --directory DIRECTORY
                        directory of structures files
  -f FILE_NAME, --file_name FILE_NAME
                        Name for saving generated features.
  -n N_JOBS, --n_jobs N_JOBS
                        Number of cpu cores for parallelization
```

Example for OIC:

```
python ./reindeer_software.py -m OIC -d ../test/ -f feature_vector_oic.csv -n -1
```
## As Module
REINDEER can also be used as a module.

Example for OIC:

```
from feature_generators import oic_dwic
from script import utils

oic = oic_dwic.InterAtomicContact(
    pathfiles="../test/",
    filename="oic_fv.csv",
    ligand_format="mol2",
    amino_acid_classes=utils.amino_acid_classes_OIC,
    cutoff=12.0,
    feature_type="OIC",
    exp=None,
)
```

## Case Study
[CaseStudy.ipynb](https://github.com/miladrayka/reindeer_software/blob/main/CaseStudy.ipynb) contains all code to reproduce the case study section of the paper on Google COLAB.

## System Specification

REINDEER is tested on the following system:

| OS  |  RAM | CPU  | Browser |
| ------------ | ------------ | ------------ |------------ |
| Windows 10  | 8.00 GB  |  AMD FX-770K Quad Core Processor (3.5 GHz) | 

We don't assume using macOS or Linux. 

## Development

To ensure code quality and consistency the following extensions of VSCode are used during development:

- black
- isort
- pylance
- pylint
- flake8
- AI python docstring generators

## Copy Right
Copyright (c) 2024, Milad Rayka
