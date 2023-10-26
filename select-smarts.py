"""
/******************************************************************************
  This source file is part of the Avogadro project.

  This source code is released under the New BSD License, (the "License").
******************************************************************************/
"""

import argparse
import json
import sys

from rdkit import Chem
from rdkit.Chem import AllChem

# Some globals:
debug = True


def getOptions():
    userOptions = {
        'SMARTS': {
            'type': 'string',
            'label': 'SMARTS Pattern',
            'default': 'a'
        }
    }

    opts = {'userOptions': userOptions }
    opts['inputMoleculeFormat'] = 'sdf'

    return opts


def runCommand():
    # Read options from stdin
    stdinStr = sys.stdin.read()

    # Parse the JSON strings
    opts = json.loads(stdinStr)
    mol = Chem.MolFromMolBlock(opts['sdf'])
    smarts = Chem.MolFromSmarts(opts['SMARTS'])
    selected = []
    for match in mol.GetSubstructMatches(smarts):
        for atom in match:
            selected.append(atom)

    # Just indicate that we want to select matching atoms
    result = {
        'selectedAtoms': selected,
        'append': True
    }
    return result

if __name__ == "__main__":
    parser = argparse.ArgumentParser('RDKit Select SMARTS Matches')
    parser.add_argument('--debug', action='store_true')
    parser.add_argument('--print-options', action='store_true')
    parser.add_argument('--run-command', action='store_true')
    parser.add_argument('--display-name', action='store_true')
    parser.add_argument('--menu-path', action='store_true')
    parser.add_argument('--lang', nargs='?', default='en')
    args = vars(parser.parse_args())

    debug = args['debug']

    if args['display_name']:
        print("Select by SMARTS Matches…")
    if args['menu_path']:
        print("&Select")
    if args['print_options']:
        print(json.dumps(getOptions()))
    elif args['run_command']:
        print(json.dumps(runCommand()))
