#!/usr/bin/env python3
"""Reproduce publication figures from retained LAMMPS checkpoint CSVs only."""
import argparse
from pathlib import Path
import sys
sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "common"))
from publicationStyle import generate

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    action = parser.add_mutually_exclusive_group(required=True)
    action.add_argument("--write", action="store_true")
    action.add_argument("--check", action="store_true")
    args = parser.parse_args()
    generate(Path(__file__).resolve().parents[1], "LAMMPS", check=args.check)
