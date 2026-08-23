#  This source file is part of the Avogadro project
#  This source code is released under the 3-Clause BSD License, (see "LICENSE").
#  https://github.com/ghutchis/avogadro-rdkit/

"""3D conformer generation and force field optimization using RDKit."""

from io import StringIO

from rdkit import Chem
from rdkit.Chem import AllChem

from .command import report_progress

# RDKit force field energies are in kcal/mol; Avogadro works in kJ/mol.
KCAL_TO_KJ = 4.184

# Matches the default of MMFFOptimizeMoleculeConfs() / UFFOptimizeMoleculeConfs().
MAX_FF_ITERATIONS = 200


class _Progress:
    """Report a sequence of steps to Avogadro as one determinate progress bar.

    Conformer generation runs in two phases (embedding, then force field), so
    the steps of both are counted against a single total to keep the bar
    monotonic instead of restarting it halfway through.
    """

    def __init__(self, total: int):
        self.total = total
        self.done = 0

    def step(self, message: str):
        report_progress(message, self.done, self.total)
        self.done += 1

    def finish(self, message: str):
        report_progress(message, self.total, self.total)


def _embed_conformers(m, n_confs: int, progress: _Progress) -> list:
    """Embed `n_confs` ETKDGv3 conformers, reporting progress for each.

    The conformers are embedded one at a time rather than in a single
    EmbedMultipleConfs() call so that progress can be reported as they appear.
    Each call keeps the conformers already on the molecule (clearConfs) and
    draws a fresh random seed (randomSeed of -1, the ETKDGv3 default), so the
    result is the same ensemble a single call would produce.
    """
    params = AllChem.ETKDGv3()
    params.clearConfs = False

    cids = []
    for i in range(n_confs):
        progress.step(f"Generating conformer {i + 1} of {n_confs}")
        cids.extend(AllChem.EmbedMultipleConfs(m, numConfs=1, params=params))
    return cids


def _conformer_energies(m, cids, ff: str, progress: _Progress) -> dict:
    """Return {confId: energy (kcal/mol)} for the requested force field.

    "MMFF94" optimizes with MMFF94 (falling back to UFF for atoms without MMFF
    parameters), "UFF" optimizes with UFF, and "None" leaves the ETKDG geometry
    untouched but still reports UFF single-point energies so the conformers can
    be ranked. Any conformer whose energy cannot be computed is omitted.
    """
    optimize = ff != "None"
    use_mmff = ff == "MMFF94" and AllChem.MMFFHasAllMoleculeParams(m)
    mmff_props = AllChem.MMFFGetMoleculeProperties(m) if use_mmff else None

    label = "Optimizing" if optimize else "Evaluating"
    total = len(cids)

    energies = {}
    for i, cid in enumerate(cids):
        progress.step(f"{label} conformer {i + 1} of {total}")
        if mmff_props is not None:
            field = AllChem.MMFFGetMoleculeForceField(m, mmff_props, confId=cid)
        else:
            field = AllChem.UFFGetMoleculeForceField(m, confId=cid)
        if field is None:
            continue
        if optimize:
            field.Minimize(maxIts=MAX_FF_ITERATIONS)
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

    # One embedding step and one force field step per conformer.
    progress = _Progress(2 * n_confs)

    cids = _embed_conformers(m, n_confs, progress)
    if not cids:
        return {}

    energies = _conformer_energies(m, cids, ff, progress)

    # Keep only conformers with a known energy so the SDF energy field stays
    # parallel to the coordinate sets, and order them lowest energy first.
    ordered = sorted((c for c in cids if c in energies),
                     key=lambda c: energies[c])
    if not ordered:
        # No energies available (e.g. UFF also lacks parameters) - still return
        # the conformers, just without energy tags.
        ordered = cids

    progress.finish(f"Writing {len(ordered)} conformers")

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
    ff = options.get("ff", "MMFF94")

    # Without a force field the conformers are not ranked, so only the
    # embedding phase has steps to report.
    progress = _Progress(n_confs if ff == "None" else 2 * n_confs)

    cids = _embed_conformers(m, n_confs, progress)
    if not cids:
        return {}

    best_cid = cids[0]
    if ff != "None":
        energies = _conformer_energies(m, cids, ff, progress)
        if energies:
            best_cid = min(energies, key=lambda c: energies[c])

    progress.finish("Building structure")

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
    report_progress("Optimizing geometry with MMFF94")
    AllChem.MMFFOptimizeMolecule(m)
    return {"moleculeFormat": "sdf", "sdf": Chem.MolToMolBlock(m)}


def uff(avo_input: dict) -> dict:
    """Optimize geometry with UFF force field."""
    m = Chem.MolFromMolBlock(avo_input["sdf"])
    m = Chem.AddHs(m)
    report_progress("Optimizing geometry with UFF")
    AllChem.UFFOptimizeMolecule(m)
    return {"moleculeFormat": "sdf", "sdf": Chem.MolToMolBlock(m)}
