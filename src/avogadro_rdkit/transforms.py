#  This source file is part of the Avogadro project
#  This source code is released under the 3-Clause BSD License, (see "LICENSE").
#  https://github.com/ghutchis/avogadro-rdkit/

"""Molecular transformations: canonicalization, tautomers, standardization, stereochemistry."""

import rdkit.Chem.rdMolTransforms
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.MolStandardize import Standardizer, rdMolStandardize


def canon_conf(avo_input: dict) -> dict:
    """Place molecule in a canonical (standard) frame of reference."""
    m = Chem.MolFromMolBlock(avo_input["sdf"])
    rdkit.Chem.rdMolTransforms.CanonicalizeMol(m)
    m = Chem.AddHs(m)
    return {"moleculeFormat": "sdf", "sdf": Chem.MolToMolBlock(m)}


def tautomer(avo_input: dict) -> dict:
    """Convert to canonical tautomer and generate a new 3D conformer."""
    m = Chem.MolFromMolBlock(avo_input["sdf"])
    enumerator = rdMolStandardize.TautomerEnumerator()
    canon = enumerator.Canonicalize(m)
    m = Chem.AddHs(canon)
    AllChem.EmbedMolecule(m, AllChem.ETKDGv3())
    return {"moleculeFormat": "sdf", "sdf": Chem.MolToMolBlock(m)}


def standardize(avo_input: dict) -> dict:
    """Standardize molecule and generate a new 3D conformer."""
    m = Chem.MolFromMolBlock(avo_input["sdf"])
    s = Standardizer()
    canon = s.standardize(m)
    m = Chem.AddHs(canon)
    AllChem.EmbedMolecule(m, AllChem.ETKDGv3())
    return {"moleculeFormat": "sdf", "sdf": Chem.MolToMolBlock(m)}


def assign_stereo(avo_input: dict) -> dict:
    """Assign stereochemistry labels (R/S) to chiral centers in the CJSON."""
    cjson = avo_input["cjson"]
    mol = Chem.MolFromMolBlock(avo_input["sdf"])

    atom_labels = cjson["atoms"].get("labels", [])
    if len(atom_labels) < mol.GetNumAtoms():
        atom_labels = [""] * mol.GetNumAtoms()

    centers = Chem.FindMolChiralCenters(mol)
    for atom_id, label in centers:
        atom_labels[atom_id] = f"({label})"

    cjson["atoms"]["labels"] = atom_labels
    return {"moleculeFormat": "cjson", "cjson": cjson}
