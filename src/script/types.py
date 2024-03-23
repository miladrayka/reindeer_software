"""Type hints for function and class annotations
"""

from typing_extensions import TypeAlias
from rdkit import Chem

Mol: TypeAlias = Chem.rdchem.Mol
