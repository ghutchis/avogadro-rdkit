#  This source file is part of the Avogadro project
#  This source code is released under the 3-Clause BSD License, (see "LICENSE").
#  https://github.com/ghutchis/avogadro-rdkit/

"""Dispatch feature identifier to the appropriate handler."""

import json
import sys


def run(args):
    avo_input = json.load(sys.stdin)
    output = None

    match args.feature:
        case "etkdg":
            from .conformers import etkdg
            output = etkdg(avo_input)
        case "mmff94":
            from .conformers import mmff94
            output = mmff94(avo_input)
        case "uff":
            from .conformers import uff
            output = uff(avo_input)
        case "canon-conf":
            from .transforms import canon_conf
            output = canon_conf(avo_input)
        case "tautomer":
            from .transforms import tautomer
            output = tautomer(avo_input)
        case "standardize":
            from .transforms import standardize
            output = standardize(avo_input)
        case "assign-stereo":
            from .transforms import assign_stereo
            output = assign_stereo(avo_input)
        case "select-smarts":
            from .selection import select_smarts
            output = select_smarts(avo_input)

    if output is not None:
        print(json.dumps(output))
