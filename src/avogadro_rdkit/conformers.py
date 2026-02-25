#  This source file is part of the Avogadro project
#  This source code is released under the 3-Clause BSD License, (see "LICENSE").
#  https://github.com/ghutchis/avogadro-rdkit/

"""3D conformer generation and force field optimization using RDKit."""

from rdkit import Chem
from rdkit.Chem import AllChem


def etkdg(avo_input: dict) -> dict:
    """Generate a 3D conformer using ETKDGv3 with optional FF optimization."""
    m = Chem.MolFromMolBlock(avo_input["sdf"])
    m = Chem.AddHs(m)
    AllChem.EmbedMolecule(m, AllChem.ETKDGv3())

    ff = avo_input.get("options", {}).get("ff", "None")
    if ff == "UFF":
        AllChem.UFFOptimizeMolecule(m)
    elif ff == "MMFF94":
        AllChem.MMFFOptimizeMolecule(m)

    return {"moleculeFormat": "sdf", "sdf": Chem.MolToMolBlock(m)}


def mmff94(avo_input: dict) -> dict:
    """Optimize geometry with MMFF94 force field."""
    m = Chem.MolFromMolBlock(avo_input["sdf"])
    m = Chem.AddHs(m)
    AllChem.MMFFOptimizeMolecule(m)
    return {"moleculeFormat": "sdf", "sdf": Chem.MolToMolBlock(m)}


def uff(avo_input: dict) -> dict:
    """Optimize geometry with UFF force field."""
    m = Chem.MolFromMolBlock(avo_input["sdf"])
    m = Chem.AddHs(m)
    AllChem.UFFOptimizeMolecule(m)
    return {"moleculeFormat": "sdf", "sdf": Chem.MolToMolBlock(m)}
