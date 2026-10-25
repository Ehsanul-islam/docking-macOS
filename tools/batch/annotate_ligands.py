#!/usr/bin/env python3
"""Add the reference ligand (largest drug-like hetero group, by the app's own classification) to each manifest item.

  python tools/batch/annotate_ligands.py --manifest M.json --data DIR
Entries with no ligand get "ligand": null, which the batch treats as 'no site defined'.
"""
import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "chemhelper"))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--manifest", required=True)
    ap.add_argument("--data", required=True)
    a = ap.parse_args()
    import gemmi

    from dockchem import hetatm
    m = json.load(open(a.manifest))
    n_lig = 0
    for it in m["items"]:
        f = next(iter(Path(a.data).glob(it["id"] + ".*")), None)
        it["ligand"], it["file"], it["n_atoms"] = None, f.name if f else None, None
        if f is None:
            continue
        try:
            st = gemmi.read_structure(str(f))
        except Exception as e:  # noqa: BLE001
            it["read_error"] = str(e)[:100]
            continue
        best = (0, None)
        for ch in st[0]:
            for r in ch:
                if r.het_flag != "H":
                    continue
                heavy = sum(1 for x in r if x.element.name != "H")
                if hetatm.classify(r.name, heavy) == "ligand" and heavy > best[0]:
                    best = (heavy, r.name)
        it["ligand"] = best[1]
        it["n_atoms"] = sum(len(r) for ch in st[0] for r in ch)
        n_lig += best[1] is not None
    json.dump(m, open(a.manifest, "w"), indent=1)
    print(f"{n_lig}/{len(m['items'])} entries have a reference ligand")


if __name__ == "__main__":
    main()
