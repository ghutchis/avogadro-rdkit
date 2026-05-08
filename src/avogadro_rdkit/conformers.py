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


def insert_smiles(avo_input: dict) -> dict:
    """Build a 3D structure from SMILES via ETKDG, returning the lowest-energy conformer."""
    options = avo_input.get("options", {})
    smiles = options.get("SMILES", "").strip()
    if not smiles:
        return {}

    m = Chem.MolFromSmiles(smiles)
    if m is None:
        return {}
    m = Chem.AddHs(m)

    n_confs = int(options.get("conformers", 25))
    cids = AllChem.EmbedMultipleConfs(m, numConfs=n_confs, params=AllChem.ETKDGv3())
    if not cids:
        return {}

    ff = options.get("ff", "MMFF94")
    if ff == "MMFF94":
        if AllChem.MMFFHasAllMoleculeParams(m):
            results = AllChem.MMFFOptimizeMoleculeConfs(m, mmffVariant="MMFF94")
        else:
            results = AllChem.UFFOptimizeMoleculeConfs(m)
    elif ff == "UFF":
        results = AllChem.UFFOptimizeMoleculeConfs(m)
    else:
        results = [(0, 0.0) for _ in cids]

    best_cid = cids[min(range(len(cids)), key=lambda i: results[i][1])]
    for cid in list(cids):
        if cid != best_cid:
            m.RemoveConformer(cid)

    return {
        "moleculeFormat": "sdf",
        "sdf": Chem.MolToMolBlock(m),
        "append": True,
    }


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
