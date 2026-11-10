#!/usr/bin/env python3
"""Reproduce publication figures from retained PFC analytical checkpoints only."""
import argparse
from pathlib import Path
import sys
sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "common"))
from publicationStyle import generate

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parent)
    action = parser.add_mutually_exclusive_group()
    action.add_argument("--write", action="store_true")
    action.add_argument("--check", action="store_true")
    args = parser.parse_args()
    generate(args.root, "PFC", check=args.check)
