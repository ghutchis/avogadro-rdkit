#  This source file is part of the Avogadro project
#  This source code is released under the 3-Clause BSD License, (see "LICENSE").
#  https://github.com/ghutchis/avogadro-rdkit/

"""Atom selection using SMARTS pattern matching."""

from rdkit import Chem


def select_smarts(avo_input: dict) -> dict:
    """Select atoms matching a SMARTS pattern."""
    mol = Chem.MolFromMolBlock(avo_input["sdf"], removeHs=False, sanitize=False)
    smarts_string = avo_input.get("options", {}).get("SMARTS", "a")
    smarts = Chem.MolFromSmarts(smarts_string)

    selected = []
    for match in mol.GetSubstructMatches(smarts):
        for atom in match:
            selected.append(atom)

    return {"selectedAtoms": selected, "append": True}
