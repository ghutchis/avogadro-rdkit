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

# Some globals:
debug = True

def getOptions():
    opts = {'userOptions': {} }
    opts['inputMoleculeFormat'] = 'sdf'

    return opts


def generate(opts):
    cjson = opts['cjson']

    # get the labels from RDKit    
    mol = Chem.MolFromMolBlock(opts['sdf'])

    # this might be empty, so expand as needed
    atomLabels = cjson['atoms']['labels']
    if len(atomLabels) < mol.GetNumAtoms():
        atomLabels = list(str('') * mol.GetNumAtoms())

    centers = Chem.FindMolChiralCenters(mol)
    # add them to the CJSON
    for id, label in centers:
        atomLabels[id] = f"({label})"

    return cjson


def runCommand():
    # Read options from stdin
    stdinStr = sys.stdin.read()

    # Parse the JSON strings
    opts = json.loads(stdinStr)

    # Replace this molecule with the new labels
    result = {}
    result['moleculeFormat'] = 'cjson'
    result['cjson'] = generate(opts)
    return result

if __name__ == "__main__":
    parser = argparse.ArgumentParser('RDKit Assign Stereo')
    parser.add_argument('--debug', action='store_true')
    parser.add_argument('--print-options', action='store_true')
    parser.add_argument('--run-command', action='store_true')
    parser.add_argument('--display-name', action='store_true')
    parser.add_argument('--menu-path', action='store_true')
    parser.add_argument('--lang', nargs='?', default='en')
    args = vars(parser.parse_args())

    debug = args['debug']

    if args['display_name']:
        print("Assign Stereo Labels")
    if args['menu_path']:
        print("&Extensions|RDKit")
    if args['print_options']:
        print(json.dumps(getOptions()))
    elif args['run_command']:
        print(json.dumps(runCommand()))
