#!/usr/bin/env python3
"""Run every wwPDB Chemical Component Dictionary ligand through ligand preparation and audit the result.

  python tools/batch/ligand_batch.py --out /tmp/ccd_run --sample 500          # quick look
  python tools/batch/ligand_batch.py --out /tmp/ccd_run                       # everything, resumable
  python tools/batch/ligand_batch.py --out /tmp/ccd_run --report              # summarise only

A ligand counts as OK only if preparation succeeded AND the output passes the audit (finite coordinates,
no overlapping atoms, PDBQT with a torsion tree, heavy-atom count consistent with the input). "Did not
crash" is not the bar.
"""
import argparse
import collections
import gzip
import json
import math
import random
import re
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "chemhelper"))

CCD = Path(__file__).resolve().parents[2] / "chemhelper" / "dockchem" / "data" / "ccd_ligands.json.gz"


def work(item):
    """Prepare one dictionary ligand (runs in a worker process)."""
    from rdkit import Chem, RDLogger
    RDLogger.DisableLog("rdApp.*")
    from dockchem.ligprep import LigandPrepOptions, prepare_ligand
    smis = item["smiles"]
    mol, used = None, None
    for k, s in enumerate(smis):
        m = Chem.MolFromSmiles(s)
        if m is not None:
            mol, used = m, k
            break
    if mol is None:
        return {"category": "dictionary_smiles_unparseable", "error": "no dictionary SMILES parses in RDKit"}
    heavy_in = mol.GetNumHeavyAtoms()
    out = {"heavy_in": heavy_in, "smiles_index": used, "type": item.get("type", "")}
    if heavy_in < 2:
        out["category"] = "skipped_atom_or_ion"
        return out
    res = prepare_ligand(smis[used], item["id"], LigandPrepOptions())
    out["variants"] = len(res)
    problems = []
    for r in res:
        if r.status == "failed":
            kind = "refused_unsupported_element" if r.error.startswith("contains ") else "prep_failed"
            problems.append((kind, r.error[:200]))
            continue
        problems += audit(r, heavy_in)
    if problems:
        out["category"] = problems[0][0]
        out["error"] = problems[0][1]
        out["all_problems"] = sorted({p[0] for p in problems})
    else:
        out["category"] = "ok"
    out["warnings"] = sorted({w[:80] for r in res for w in r.warnings})[:5]
    return out


def audit(r, heavy_in):
    from rdkit import Chem
    problems = []
    mol = r.mol
    conf = mol.GetConformer()
    xyz = [conf.GetAtomPosition(i) for i in range(mol.GetNumAtoms())]
    if any(not all(math.isfinite(c) for c in (p.x, p.y, p.z)) for p in xyz):
        return [("bad_coordinates", "non-finite coordinates")]
    dmin = min((xyz[i].Distance(xyz[j]) for i in range(len(xyz)) for j in range(i + 1, len(xyz))), default=9)
    if dmin < 0.6:
        problems.append(("overlapping_atoms", f"two atoms {dmin:.2f} A apart"))
    pt = Chem.GetPeriodicTable()
    bad_bonds = [(b.GetBeginAtom().GetSymbol(), b.GetEndAtom().GetSymbol(),
                  round(xyz[b.GetBeginAtomIdx()].Distance(xyz[b.GetEndAtomIdx()]), 2))
                 for b in mol.GetBonds()
                 if not 0.85 <= xyz[b.GetBeginAtomIdx()].Distance(xyz[b.GetEndAtomIdx()]) <= 1.5 * (
                     pt.GetRcovalent(b.GetBeginAtom().GetAtomicNum()) + pt.GetRcovalent(b.GetEndAtom().GetAtomicNum()))]
    if bad_bonds:
        problems.append(("implausible_bond_length", f"{len(bad_bonds)} bond(s) implausible, e.g. {bad_bonds[0]}"))
    heavy_out = sum(1 for a in mol.GetAtoms() if a.GetAtomicNum() > 1)
    stripped = any("salt" in s.lower() or "fragment" in s.lower() for s in r.steps + r.warnings)
    if heavy_out != heavy_in and not stripped:
        problems.append(("heavy_atom_count_changed", f"{heavy_in} -> {heavy_out}"))
    if not r.pdbqt:
        problems.append(("no_pdbqt", "empty PDBQT"))
    else:
        atoms = [l for l in r.pdbqt.splitlines() if l.startswith(("ATOM", "HETATM"))]
        if not atoms or "ROOT" not in r.pdbqt or "TORSDOF" not in r.pdbqt:
            problems.append(("malformed_pdbqt", "missing ATOM/ROOT/TORSDOF"))
        elif len(atoms) < heavy_out:
            problems.append(("pdbqt_atoms_missing", f"{len(atoms)} PDBQT atoms for {heavy_out} heavy atoms"))
        else:
            bad = [l for l in atoms if not all(math.isfinite(float(l[30 + 8 * k:38 + 8 * k])) for k in range(3))]
            if bad:
                problems.append(("bad_pdbqt_coordinates", "non-finite"))
    return problems


def load_items(sample, seed, kinds, only_ids=None):
    data = json.load(gzip.open(CCD))["entries"]
    items = [{"id": k, "smiles": v["s"], "type": v.get("t", "")} for k, v in data.items() if v.get("s")]
    if kinds:
        items = [i for i in items if any(t in i["type"].upper() for t in kinds)]
    if only_ids is not None:
        items = [i for i in items if i["id"] in only_ids]
    items.sort(key=lambda i: i["id"])
    if sample:
        random.Random(seed).shuffle(items)
        items = items[:sample]
    return items


def norm_error(e):
    return re.sub(r"\d+", "N", e or "")[:70]


def report(path):
    rows = [json.loads(l) for l in Path(path).read_text().splitlines() if l.strip()]
    cat = collections.Counter(r["category"] for r in rows)
    n = len(rows)
    print(f"\n{n} ligands")
    for c, k in cat.most_common():
        print(f"  {c:34s}{k:7d}  {100 * k / n:5.1f}%")
    scored = [r for r in rows if r["category"] not in ("skipped_atom_or_ion", "dictionary_smiles_unparseable", "refused_unsupported_element")]
    ok = sum(r["category"] == "ok" for r in scored)
    print(f"\nprepared and audited OK: {ok}/{len(scored)} = {100 * ok / max(1, len(scored)):.2f}% of ligands that are in scope (organic elements, parseable)")
    t = sorted(r.get("seconds", 0) for r in rows)
    if t:
        print(f"seconds per ligand: median {t[len(t) // 2]:.2f}, p95 {t[int(len(t) * .95)]:.2f}, max {t[-1]:.1f}")
    errs = collections.defaultdict(list)
    for r in rows:
        if r["category"] not in ("ok", "skipped_atom_or_ion"):
            errs[(r["category"], norm_error(r.get("error")))].append(r["id"])
    print("\nfailure groups (category / message): count, examples")
    for (c, e), ids in sorted(errs.items(), key=lambda kv: -len(kv[1]))[:25]:
        print(f"  {len(ids):5d}  {c} / {e}   e.g. {', '.join(ids[:4])}")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("--sample", type=int, default=0)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--timeout", type=float, default=120)
    ap.add_argument("--kinds", nargs="*", help="only CCD types containing these words, e.g. NON-POLYMER")
    ap.add_argument("--report", action="store_true")
    ap.add_argument("--rerun-failures-from", help="JSONL of an earlier run: only ligands that were not ok there")
    a = ap.parse_args()
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    jl = out / "ligands.jsonl"
    only = None
    if a.rerun_failures_from:
        only = {r["id"] for r in map(json.loads, open(a.rerun_failures_from)) if r["category"] not in ("ok", "skipped_atom_or_ion", "dictionary_smiles_unparseable")}
    if not a.report:
        from runner import run_batch
        print(run_batch(load_items(a.sample, a.seed, a.kinds, only), "ligand_batch:work", jl, a.workers, a.timeout))
    report(jl)
