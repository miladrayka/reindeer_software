# REINDEER Software
REINDEER is a software for the sturcture-based protein-ligand feature generation.
Currently, REINDEER provides only four below feature vector:

1- Occurrence of Interatomic Contact (OIC) - **[Ref](https://academic.oup.com/bioinformatics/article/26/9/1169/199938?login=false)**

2- Distance-Weighted Interatomic Contact (DWIC) - **[Ref](https://onlinelibrary.wiley.com/doi/abs/10.1002/minf.202060084)**

3- Extended Connectivity Interaction Feature (ECIF) - **[Ref](https://academic.oup.com/bioinformatics/article/37/10/1376/5998664?login=false)**

4- Multi-Shell Occurrence of Interatomic Contact (MS-OIC) - **[Ref](https://www.frontiersin.org/articles/10.3389/fchem.2021.753002/full)**


![Logo](https://github.com/miladrayka/reindeer_software/blob/main/reindeer/logo/Logo.png)

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
REINDEER provides a GUI to make for feature generation methods.
---
![pic01](https://github.com/miladrayka/reindeer_software/blob/main/images/pic01.PNG)

---
![pic02](https://github.com/miladrayka/reindeer_software/blob/main/images/pic02.PNG)

---
![pic03](https://github.com/miladrayka/reindeer_software/blob/main/images/pic03.PNG)

---
![pic04](https://github.com/miladrayka/reindeer_software/blob/main/images/pic04.PNG)

---
![pic05](https://github.com/miladrayka/reindeer_software/blob/main/images/pic05.PNG)

---
![pic06](https://github.com/miladrayka/reindeer_software/blob/main/images/pic06.PNG)

