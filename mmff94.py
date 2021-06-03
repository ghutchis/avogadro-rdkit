#  This source file is part of the Avogadro project
#  This source code is released under the 3-Clause BSD License, (see "LICENSE").
#  https://github.com/ghutchis/avogadro-rdkit/

import argparse
import json
import sys

from rdkit import Chem
from rdkit.Chem import AllChem

# Some globals:
debug = True


def getOptions():
    opts = {}
    opts['inputMoleculeFormat'] = 'sdf'

    return opts


def generate(opts):
    m = Chem.MolFromMolBlock(opts['sdf'])
    m = Chem.AddHs(m)

    AllChem.MMFFOptimizeMolecule(m)

    return Chem.MolToMolBlock(m)


def runCommand():
    # Read options from stdin
    stdinStr = sys.stdin.read()

    # Parse the JSON strings
    opts = json.loads(stdinStr)

    # Replace this molecule with a new conformer in SDF
    result = {}
    result['moleculeFormat'] = 'sdf'
    result['sdf'] = generate(opts)
    return result

if __name__ == "__main__":
    parser = argparse.ArgumentParser('RDKit MMFF')
    parser.add_argument('--debug', action='store_true')
    parser.add_argument('--print-options', action='store_true')
    parser.add_argument('--run-command', action='store_true')
    parser.add_argument('--display-name', action='store_true')
    parser.add_argument('--menu-path', action='store_true')
    parser.add_argument('--lang', nargs='?', default='en')
    args = vars(parser.parse_args())

    debug = args['debug']

    if args['display_name']:
        print("MMFF94 Optimize")
    if args['menu_path']:
        print("&Extensions|RDKit")
    if args['print_options']:
        print(json.dumps(getOptions()))
    elif args['run_command']:
        print(json.dumps(runCommand()))
