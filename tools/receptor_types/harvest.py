#!/usr/bin/env python3
"""Build chemhelper/dockchem/data/receptor_atom_types.json from corpora of Meeko receptor preparations.

Input : directories each holding receptor_prepared.pdb + receptor.pdbqt (a prep run's output).
Output: {"types": {"RES": {"ATOM": "AD4type"}}, "metals": ..., "provenance": {...}}

Atoms present in the PDBQT are typed; hydrogens absent from it were merged by Meeko (non-polar). Disagreement
between two corpus entries for the same (residue, atom) is reported and the key is left out, except the
known context rule for CYX SG (disulfide -> S, free -> SA) which the writer applies itself.

  python tools/receptor_types/harvest.py DIR [DIR ...] [--out PATH]
"""
import argparse
import collections
import json
import sys
from pathlib import Path

OUT = Path(__file__).resolve().parents[2] / "chemhelper" / "dockchem" / "data" / "receptor_atom_types.json"


def pairs(dirs):
    for d in dirs:
        for pdb in Path(d).rglob("receptor_prepared.pdb"):
            q = pdb.with_name("receptor.pdbqt")
            if q.exists():
                yield pdb, q


def key(l):
    return (l[17:20].strip(), l[12:16].strip())


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("dirs", nargs="+")
    ap.add_argument("--out", default=str(OUT))
    a = ap.parse_args()
    seen = collections.defaultdict(collections.Counter)
    merged = collections.defaultdict(int)
    n = 0
    for pdb, q in pairs(a.dirs):
        n += 1
        typed = set()
        for l in q.read_text().splitlines():
            if l.startswith(("ATOM", "HETATM")):
                seen[key(l)][l[77:79].strip()] += 1
                typed.add((l[21], l[22:27].strip(), l[12:16].strip()))
        for l in pdb.read_text().splitlines():
            if l.startswith(("ATOM", "HETATM")) and l[76:78].strip() == "H" and \
                    (l[21], l[22:27].strip(), l[12:16].strip()) not in typed:
                merged[key(l)] += 1
    types, conflicts = collections.defaultdict(dict), {}
    for (res, atom), c in seen.items():
        if len(c) == 1:
            types[res][atom] = next(iter(c))
        elif (res, atom) == ("CYX", "SG"):
            types[res][atom] = "SA"          # writer applies S when the SG is in a disulfide
        else:
            conflicts[f"{res}.{atom}"] = dict(c)
    merged_h = collections.defaultdict(list)
    for (res, atom), k in merged.items():
        if (res, atom) not in seen:
            merged_h[res].append(atom)
    out = {"types": dict(types), "merged_hydrogens": {r: sorted(v) for r, v in merged_h.items()},
           "provenance": {"source": "Meeko receptor PDBQT output", "corpus_pairs": n,
                          "atoms_observed": sum(sum(c.values()) for c in seen.values())}}
    Path(a.out).write_text(json.dumps(out, indent=1, sort_keys=True))
    print(f"{n} pairs, {len(seen)} (residue, atom) keys, {len(types)} residue names; conflicts: {conflicts or 'none'}")


if __name__ == "__main__":
    sys.exit(main())
