#  This source file is part of the Avogadro project
#  This source code is released under the 3-Clause BSD License, (see "LICENSE").
#  https://github.com/ghutchis/avogadro-rdkit/

"""Entry point for the avogadro-rdkit plugin."""

import argparse

from .run import run


def main():
    parser = argparse.ArgumentParser("avogadro-rdkit")
    parser.add_argument("feature")
    parser.add_argument("--lang", nargs="?", default="en")
    parser.add_argument("--debug", action="store_true")
    args = parser.parse_args()
    run(args)
