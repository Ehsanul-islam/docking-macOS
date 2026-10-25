#!/usr/bin/env python3
"""Distil the wwPDB Chemical Component Dictionary into a compact bundle for the app.

Input : components.cif.gz  (https://files.wwpdb.org/pub/pdb/data/monomers/components.cif.gz, public domain)
Output: chemhelper/dockchem/data/ccd_ligands.json.gz  -> {id: {n: name, f: formula, t: type, s: [smiles...]}}

Only what bond-order assignment needs is kept (flat SMILES, no coordinates). SMILES preference: OpenEye
"SMILES" (flat), then ACDLabs "SMILES" (explicit double bonds); several are kept as fallbacks because RDKit
cannot parse every one (organometallics in particular).

  python tools/ligand_dictionary/build_ccd.py --ccd /tmp/components.cif.gz --last-modified "Sat, 17 Oct 2026 ..."
"""
import argparse
import gzip
import json
import sys
import time
from pathlib import Path

import gemmi

PREFERENCE = [("SMILES", "OpenEye OEToolkits"), ("SMILES", "ACDLabs"), ("SMILES_CANONICAL", "OpenEye OEToolkits"),
              ("SMILES", "CACTVS")]


def blocks(path):
    """Yield the text of each data_ block without loading the whole file."""
    cur = []
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as f:
        for line in f:
            if line.startswith("data_") and cur:
                yield "".join(cur)
                cur = []
            cur.append(line)
    if cur:
        yield "".join(cur)


def distil(text: str):
    doc = gemmi.cif.read_string(text)
    b = doc.sole_block()
    cid = b.name
    pair = lambda tag: (b.find_value(tag) or "").strip('"\' ')
    smi = {}
    t = b.find(["_pdbx_chem_comp_descriptor.type", "_pdbx_chem_comp_descriptor.program",
                "_pdbx_chem_comp_descriptor.descriptor"])
    for row in t:
        smi.setdefault((gemmi.cif.as_string(row[0]), gemmi.cif.as_string(row[1])), gemmi.cif.as_string(row[2]))
    chosen = [smi[k] for k in PREFERENCE if k in smi and smi[k]]
    if not chosen:
        return cid, None
    seen, uniq = set(), []
    for s in chosen:
        if s not in seen:
            seen.add(s); uniq.append(s)
    return cid, {"n": pair("_chem_comp.name"), "f": pair("_chem_comp.formula"), "t": pair("_chem_comp.type"),
                 "r": pair("_chem_comp.pdbx_release_status"), "s": uniq}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ccd", required=True)
    ap.add_argument("--out", default=str(Path(__file__).resolve().parents[2] / "chemhelper/dockchem/data/ccd_ligands.json.gz"))
    ap.add_argument("--last-modified", default="")
    a = ap.parse_args()
    entries, skipped, t0 = {}, 0, time.time()
    for text in blocks(a.ccd):
        try:
            cid, e = distil(text)
        except Exception:  # noqa: BLE001 - one bad block must not lose the other 40k
            skipped += 1
            continue
        if e:
            entries[cid] = e
        else:
            skipped += 1
    doc = {"source": "wwPDB Chemical Component Dictionary (components.cif.gz), public domain",
           "url": "https://files.wwpdb.org/pub/pdb/data/monomers/components.cif.gz",
           "last_modified": a.last_modified, "built": time.strftime("%Y-%m-%d"),
           "count": len(entries), "entries": entries}
    with gzip.open(a.out, "wt", encoding="utf-8", compresslevel=9) as f:
        json.dump(doc, f, separators=(",", ":"), sort_keys=True)
    print(f"{len(entries)} components, {skipped} without SMILES/unparseable, {time.time() - t0:.0f}s, "
          f"{Path(a.out).stat().st_size / 1e6:.1f} MB -> {a.out}")


if __name__ == "__main__":
    main()
