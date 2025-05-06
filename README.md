[![Python 3.9](https://img.shields.io/badge/python-3.9-blue.svg)](https://www.python.org/downloads/release/python-390/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
![Last Commit](https://img.shields.io/github/last-commit/miladrayka/reindeer_software)
![Issues](https://img.shields.io/github/issues/miladrayka/reindeer_software)

# 🦌 REINDEER: A Protein-Ligand Feature Generator

Efficiently generate structure-based protein-ligand interaction features for machine learning models.

![Logo](https://github.com/miladrayka/reindeer_software/blob/main/reindeer/logo/Logo.png)

---
## 📚 Table of Contents
- [Features](#-features)
- [Installation](#-installation)
- [Notes](#️-notes)
- [Usage](#️-usage)
  - [GUI](#-gui)
  - [CLI](#-cli)
  - [Python API](#-python-api)
- [Citation](#-citation)
- [Case Study](#-case-study)
- [System Specification](#-system-specification)
- [Development](#-development)
- [Original Repositories](#-original-repositories)
- [License](#-license)
- [Contact](#-contact)

---
## 🔬 Features

REINDEER supports the following feature generation methods:

- **OIC**: Occurrence of Interatomic Contact ([Ref](https://academic.oup.com/bioinformatics/article/26/9/1169/199938?login=false))
- **DWIC**: Distance-Weighted Interatomic Contact ([Ref](https://onlinelibrary.wiley.com/doi/abs/10.1002/minf.202060084))
- **ECIF**: Extended Connectivity Interaction Features ([Ref](https://academic.oup.com/bioinformatics/article/37/10/1376/5998664?login=false))
- **MS-OIC**: Multi-Shell Occurrence of Interatomic Contact ([Ref](https://www.frontiersin.org/articles/10.3389/fchem.2021.753002/full))

---
## 🚀 Installation

1. Create and activate a Python 3.9 virtual environment:
```bash
   python -m venv env
   source env/bin/activate  # On Windows: .\env\Scripts\activate
```

2. Clone the repository:

```bash
git clone https://github.com/miladrayka/reindeer_software.git
cd reindeer_software

```
3. Install dependencies:
```bash
pip install -r requirements.txt
```

---
💡 Notes
- Input protein-ligand complexes **must include hydrogen atoms**.
- Supported file formats:
	- Protein: `.pdb`
	- Ligand: `.mol2` (use `.sdf` for ECIF)

```bash
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
│   ├── 1a30_ligand.mol2
│   ├── 1a30_ligand.sdf
│   └── 1a30_protein.pdb
```

---
### 🔳 GUI

Launch the GUI:

```bash
python ./gui_launcher.py
```

👉 See the [Tutorial](https://github.com/miladrayka/reindeer_software/blob/main/Tutorial.pdf) for step-by-step instructions.

![GUI](https://github.com/miladrayka/reindeer_software/blob/main/GUI_img.PNG)

---
### 💻 CLI

Show help message:

```bash
python ./reindeer_software.py -h
```

Example (OIC feature generation):

```bash
python ./reindeer_software.py -m OIC -d ./test/ -f feature_vector_oic.csv -n -1
```

CLI Arguments:

```bash
usage: reindeer_software.py [-h] -m METHOD -d DIRECTORY -f FILE_NAME -n N_JOBS

Generate features for set of given structures

optional arguments:
  -h, --help            show this help message and exit
  -m METHOD, --method METHOD
                        Feature generation method. Only OIC, DWIC, ECIF, and MS-OIC are implemented for now.
  -d DIRECTORY, --directory DIRECTORY
                        Directory of structure files
  -f FILE_NAME, --file_name FILE_NAME
                        Name for saving generated features
  -n N_JOBS, --n_jobs N_JOBS
                        Number of CPU cores for parallelization
```

---
### ⚙️ Python API

Example (OIC feature generation):

```bash
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

---
## 📄 Citation

The paper is currently under review. Please cite the preprint:

> Rayka, M. _REINDEER: A Protein-Ligand Feature Generator Software for Machine Learning Algorithms._ [ChemRxiv](https://chemrxiv.org/engage/chemrxiv/article-details/6613916c21291e5d1d5cd171)

---
## :bar_chart: Case Study
All codes to reproduce the case study will be available after publication.

---
## 🖥️ System Specification

REINDEER is tested on:

|OS|RAM|CPU|
|---|---|---|
|Ubuntu 24.04|64.00 GB|Intel(R) Core(TM) i7-14700KF|

macOS and Windows are expected to work as well, though not officially tested.

---
## 🔧 Development

VSCode extensions used during development:

- `black`
    
- `isort`
    
- `pylance`
    
- `pylint`
    
- `flake8`
    
- AI Python docstring generators

---

## 📚 Original Repositories

- [ECIF Github](https://github.com/DIFACQUIM/ECIF)
    
- [ET-Score and RF-Score](https://github.com/miladrayka/ENS_Score)
    
- [OnionNet-2](https://github.com/zchwang/OnionNet-2)

---
## 📝 License

This project is licensed under the [MIT License](https://opensource.org/licenses/MIT).

---
## 📬 Contact

Milad Rayka  
📧 [milad.rayka@yahoo.com](mailto:milad.rayka@yahoo.com)
