#  This source file is part of the Avogadro project
#  This source code is released under the 3-Clause BSD License, (see "LICENSE").
#  https://github.com/ghutchis/avogadro-rdkit/

"""Molecular transformations: canonicalization, tautomers, standardization, stereochemistry."""

import rdkit.Chem.rdMolTransforms
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.MolStandardize import rdMolStandardize


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
    canon = rdMolStandardize.Normalize(m)
    m = Chem.AddHs(canon)
    AllChem.EmbedMolecule(m, AllChem.ETKDGv3())
    return {"moleculeFormat": "sdf", "sdf": Chem.MolToMolBlock(m)}


def assign_stereo(avo_input: dict) -> dict:
    """Assign stereochemistry labels (R/S) to chiral centers in the CJSON."""
    cjson = avo_input["cjson"]
    mol = Chem.MolFromMolBlock(avo_input["sdf"], removeHs=False)
    Chem.AssignStereochemistryFrom3D(mol)

    atom_labels = cjson["atoms"].get("labels", [])
    if len(atom_labels) < mol.GetNumAtoms():
        atom_labels = [""] * mol.GetNumAtoms()

    centers = Chem.FindMolChiralCenters(mol, useLegacyImplementation=False)
    for atom_id, label in centers:
        # remove parens to just show R/S
        label.replace("(", "").replace(")", "")
        atom_labels[atom_id] = f"({label})"

    cjson["atoms"]["labels"] = atom_labels

    # add E / Z for double bonds too
    bond_labels = cjson["bonds"].get("labels", [])
    if len(bond_labels) < mol.GetNumBonds():
        bond_labels = [""] * mol.GetNumBonds()

    for bond in mol.GetBonds():
        if bond.GetStereo() == Chem.BondStereo.STEREOE:
            bond_labels[bond.GetIdx()] = "E"
        elif bond.GetStereo() == Chem.BondStereo.STEREOZ:
            bond_labels[bond.GetIdx()] = "Z"

    cjson["bonds"]["labels"] = bond_labels

    return {"moleculeFormat": "cjson", "cjson": cjson}
