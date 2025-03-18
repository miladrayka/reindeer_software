[![Python 3.9](https://img.shields.io/badge/python-3.9-blue.svg)](https://www.python.org/downloads/release/python-360/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

# REINDEER Software

REINDEER is a software for structure-based protein-ligand feature generation.

![Logo](https://github.com/miladrayka/reindeer_software/blob/main/reindeer/logo/Logo.png)

Currently, REINDEER provides only four feature vectors:

1- Occurrence of Interatomic Contact (OIC) - **[Ref](https://academic.oup.com/bioinformatics/article/26/9/1169/199938?login=false)**

2- Distance-Weighted Interatomic Contact (DWIC) - **[Ref](https://onlinelibrary.wiley.com/doi/abs/10.1002/minf.202060084)**

3- Extended Connectivity Interaction Feature (ECIF) - **[Ref](https://academic.oup.com/bioinformatics/article/37/10/1376/5998664?login=false)**

4- Multi-Shell Occurrence of Interatomic Contact (MS-OIC) - **[Ref](https://www.frontiersin.org/articles/10.3389/fchem.2021.753002/full)**

## Citation

The paper is currently under review. For now, you can cite the below paper:

[REINDEER: A Protein-Ligand Feature Generator Software for Machine Learning Algorithms](https://chemrxiv.org/engage/chemrxiv/article-details/6613916c21291e5d1d5cd171)

## Contact

Milad Rayka, milad.rayka@yahoo.com

## Install

1- First install Python (3.9) then make a virtual environment and activate it.
```
python -m venv env
.\env\Scripts\activate
```
Which *env* is the location to create the virtual environment.

2- Clone *reindeer_software* Github repository.
```
git clone https://github.com/miladrayka/reindeer_software.git
```
3- Change your directory to *reindeer_software*.

4- Install required packages with pip.
```
pip install -r requirements.txt
```
## Notes

1- Provided protein-ligand complex should have hydrogen atoms

2- File formats for protein and ligand are *.pdb* and *.mol2*. 
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

## Usage
REINDEER provides GUI, CLI, and using within Python codes for feature generation.

### Graphical User Interface (GUI)
After changing your directory to *reindeer_software* type the following code for running GUI:
```
python ./gui_launcher.py
```
For example, check the [Tutorial](https://github.com/miladrayka/reindeer_software/blob/main/Tutorial.pdf) file.

![GUI](https://github.com/miladrayka/reindeer_software/blob/main/GUI_img.PNG)

### Command Line Interface (CLI)
For access to CLI, type the following command (you should be in *reindeer_software* directory):
```
python ./reindeer_software.py -h
```
The output is like this:

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
python ./reindeer_software.py -m OIC -d ./test/ -f feature_vector_oic.csv -n -1
```
### Within Python
REINDEER can also be used within Python codes.

Example for OIC:

```
from reindeer.feature_generators import oic_dwic
from reindeer.script import utils

oic = oic_dwic.InterAtomicContact(
    pathfiles="./test/",
    filename="oic_fv.csv",
    ligand_format="mol2",
    amino_acid_classes=utils.amino_acid_classes_OIC,
    cutoff=12.0,
    feature_type="OIC",
    exp=None,
)

oic.generate_features(n_jobs=-1)
```
## Case Study
All codes to reproduce the case study will be available in a Zenodo repository after publication.

## System Specification

REINDEER is tested on the following system:

| OS  |  RAM | CPU  |
| ------------ | ------------ | ------------ |
| Ubuntu 24.04 | 64.00 GB  |  Intel(R) Core(TM) i7-14700KF | 

We don't assume using macOS or Windows can cause a problem. 

## Development

To ensure code quality and consistency the following extensions of VSCode are used during development:

- black
- isort
- pylance
- pylint
- flake8
- AI Python docstring generators

## Original Repository

The following repositories were used for the development of REINDEER:

* [ECIF Github](https://github.com/DIFACQUIM/ECIF) for ECIF method.

* [ET-Score and RF-Score](https://github.com/miladrayka/ENS_Score) for DWIC and OIC methods.

* [OnionNet-2](https://github.com/zchwang/OnionNet-2) MS-OIC method.

## Copy Right
Copyright (c) 2025, Milad Rayka
