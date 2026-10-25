#!/usr/bin/env python3
"""Build a protein test manifest from the PoseBusters paper data (Zenodo 8278563): the Astex Diverse Set and the
PoseBusters benchmark set, the two docking benchmarks in widest use. Only the entry IDs and ligand codes are
taken from the archive; structures are downloaded from the PDB itself so the test inputs are what users get.

  python tools/batch/make_manifest.py --zip posebusters_paper_data.zip --n 300 --out tools/batch/sets/widely_used_300.json
"""
import argparse
import json
import zipfile


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--zip", required=True)
    ap.add_argument("--n", type=int, default=300)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    z = zipfile.ZipFile(a.zip)
    sets = {}
    for n in z.namelist():
        p = n.split("/")
        if len(p) >= 2 and p[1] and p[0] in ("astex_diverse_set", "posebusters_benchmark_set"):
            sets.setdefault(p[0], set()).add(p[1])
    seen, items = set(), []

    def add(source, entries):
        for e in sorted(entries):
            pdb, _, lig = e.partition("_")
            if pdb in seen:
                continue
            seen.add(pdb)
            items.append({"id": pdb.upper(), "ligand": lig, "source": source})

    add("astex_diverse_set", sets["astex_diverse_set"])
    pb = [i for i in sorted(sets["posebusters_benchmark_set"]) if i.split("_")[0] not in seen]
    want = max(0, a.n - len(items))
    step = len(pb) / max(1, want)
    picked = [pb[int(i * step)] for i in range(want)]          # evenly spread over the sorted IDs
    add("posebusters_benchmark_set", picked)
    json.dump({"description": "Astex Diverse Set + evenly spaced PoseBusters benchmark entries",
               "source": "https://zenodo.org/records/8278563", "items": items[:a.n]}, open(a.out, "w"), indent=1)
    print(len(items[:a.n]), "entries ->", a.out)


if __name__ == "__main__":
    main()
