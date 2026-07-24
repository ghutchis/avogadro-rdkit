#  This source file is part of the Avogadro project
#  This source code is released under the 3-Clause BSD License, (see "LICENSE").
#  https://github.com/ghutchis/avogadro-rdkit/

"""3D conformer generation and force field optimization using RDKit."""

from io import StringIO

from rdkit import Chem
from rdkit.Chem import AllChem

# RDKit force field energies are in kcal/mol; Avogadro works in kJ/mol.
KCAL_TO_KJ = 4.184


def _conformer_energies(m, cids, ff: str) -> dict:
    """Return {confId: energy (kcal/mol)} for the requested force field.

    "MMFF94" optimizes with MMFF94 (falling back to UFF for atoms without MMFF
    parameters), "UFF" optimizes with UFF, and "None" leaves the ETKDG geometry
    untouched but still reports single-point energies so the conformers can be
    ranked. Any conformer whose energy cannot be computed is omitted.
    """
    use_mmff = ff == "MMFF94" and AllChem.MMFFHasAllMoleculeParams(m)

    # Optimize in place unless the user asked for no force field.
    if ff == "MMFF94" and use_mmff:
        results = AllChem.MMFFOptimizeMoleculeConfs(m, mmffVariant="MMFF94")
        return {cid: e for cid, (_conv, e) in zip(cids, results)}
    if ff == "UFF" or (ff == "MMFF94" and not use_mmff):
        results = AllChem.UFFOptimizeMoleculeConfs(m)
        return {cid: e for cid, (_conv, e) in zip(cids, results)}

    # ff == "None": single-point energies only (MMFF94 if available, else UFF).
    energies = {}
    mmff_props = AllChem.MMFFGetMoleculeProperties(m) if use_mmff else None
    for cid in cids:
        if mmff_props is not None:
            field = AllChem.MMFFGetMoleculeForceField(m, mmff_props, confId=cid)
        else:
            field = AllChem.UFFGetMoleculeForceField(m, confId=cid)
        if field is not None:
            energies[cid] = field.CalcEnergy()
    return energies


def etkdg(avo_input: dict) -> dict:
    """Generate multiple 3D conformers with ETKDGv3 and optional FF optimization.

    Returns every conformer (lowest energy first) as a multi-conformer SDF, each
    block tagged with an ``<energy>`` data field so Avogadro can plot the
    relative conformer energies.
    """
    m = Chem.MolFromMolBlock(avo_input["sdf"])
    if m is None:
        return {}
    m = Chem.AddHs(m)

    options = avo_input.get("options", {})
    n_confs = max(1, int(options.get("conformers", 10)))
    ff = options.get("ff", "MMFF94")

    cids = list(AllChem.EmbedMultipleConfs(m, numConfs=n_confs,
                                           params=AllChem.ETKDGv3()))
    if not cids:
        return {}

    energies = _conformer_energies(m, cids, ff)

    # Keep only conformers with a known energy so the SDF energy field stays
    # parallel to the coordinate sets, and order them lowest energy first.
    ordered = sorted((c for c in cids if c in energies),
                     key=lambda c: energies[c])
    if not ordered:
        # No energies available (e.g. UFF also lacks parameters) - still return
        # the conformers, just without energy tags.
        ordered = cids

    sio = StringIO()
    writer = Chem.SDWriter(sio)
    for cid in ordered:
        if cid in energies:
            m.SetDoubleProp("energy", energies[cid] * KCAL_TO_KJ)
        elif m.HasProp("energy"):
            m.ClearProp("energy")
        writer.write(m, confId=cid)
    writer.close()

    return {"moleculeFormat": "sdf", "sdf": sio.getvalue()}


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
