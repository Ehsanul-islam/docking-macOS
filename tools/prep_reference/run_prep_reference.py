#!/usr/bin/env python3
"""Run the receptor-prep pipeline over a set of real PDB entries (plan §12.2 groundwork).

Each entry runs in its own process with a timeout, so one failure never stops the batch.
Output: <out>/<ID>/ (prepared files) plus <out>/summary.json and summary.md.

  python tools/prep_reference/run_prep_reference.py --pdb-dir pdb_structures_50 --out /tmp/prep_ref
"""
import argparse
import json
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

VARIANTS = {"HID", "HIE", "HIP", "ASH", "GLH", "LYN", "CYX", "CYM"}   # protonation-state residue names

# Entries whose INPUT is defective, so a degraded result is the correct outcome (reason documented).
EXPECTED_DEGRADED = {
    "1TUB": "input has overlapping heavy atoms (GLU B71 CD is 0.4 A from ALA B99 CB)",
}

# id -> reference ligand residue name (None = apo / no ligand of interest)
SET = {
    "1IEP": "STI", "1M17": "AQ4", "1STP": "BTN", "3PTB": "BEN", "4G5J": "0WN", "1K4C": "F09",
    "2RH1": "CAU", "1AKE": "AP5", "1ATP": "ATP", "1GPB": "PLP", "3SN6": "P0G", "1MBN": "HEM",
    "1HHO": "HEM", "1TUB": "TXL", "4MDH": "NAD", "1JFF": "TA1",
    "1CA2": None, "3CPA": None, "2PTC": None, "1CLL": None,
    "1CRN": None, "1UBQ": None, "1AKI": None, "1L2Y": None, "5RSA": None, "1EMA": None, "1TIM": None,
    "4INS": None,
}


def _residue_keys(path: str, fmt: str, chains=None) -> dict:
    """{(chain, seq+icode): resname} for amino-acid residues."""
    out = {}
    if fmt == "pdb":
        import gemmi
        st = gemmi.read_structure(path)
        for ch in st[0]:
            if chains is not None and ch.name not in chains:
                continue
            for r in ch:
                if r.het_flag == "H" and fmt == "pdb" and chains is not None:
                    continue            # bound peptides (HETATM) are ligands, not receptor residues
                info = gemmi.find_tabulated_residue(r.name)
                if r.name in VARIANTS or (info is not None and info.is_amino_acid()):
                    out[(ch.name, f"{r.seqid.num}{r.seqid.icode.strip()}")] = r.name
    else:      # pdbqt
        for l in open(path):
            if l.startswith("ATOM"):
                out[(l[21], l[22:27].strip())] = l[17:20].strip()
    return out


def worker(pdb: str, ref: str | None, outdir: str) -> None:
    from dockchem.recprep import ReceptorPrepOptions, prepare_receptor, _meeko_templates
    res = prepare_receptor(pdb, outdir, ReceptorPrepOptions(reference_ligand=ref))
    r = res.report
    src = _residue_keys(pdb, "pdb", set(r["chains"]))
    prep = _residue_keys(str(res.pdb), "pdb")
    qt = _residue_keys(str(res.pdbqt), "pdbqt")
    templated = _meeko_templates()
    lost_pdb = sorted(f"{c}:{n}{src[(c, n)]}" for (c, n) in src if (c, n) not in prep)
    # a residue missing from the PDBQT is only acceptable when Meeko has no template for its name
    lost_qt = sorted(f"{c}:{n}{prep[(c, n)]}" for (c, n) in prep
                     if (c, n) not in qt and prep[(c, n)] in templated)
    print(json.dumps({
        "atoms_before": r["atoms_before"], "atoms_after": r["atoms_after"], "chains": r["chains"],
        "n_changes": len(r["changes"]), "n_warnings": len(r["warnings"]), "warnings": r["warnings"],
        "ligands": [{k: v for k, v in l.items() if k != "path"} for l in res.ligands],
        "pdbqt_atoms": res.pdbqt.read_text().count("ATOM") + res.pdbqt.read_text().count("HETATM"),
        "stages": sorted({c["stage"] for c in r["changes"]}),
        "protonation": [c["change"] for c in r["changes"] if c["stage"] == "protonation"],
        "lost_residues_pdb": lost_pdb, "lost_residues_pdbqt": lost_qt,
        "status_detail": r["status"], "n_dropped": len(r["dropped_residues"]),
    }))


def run_one(pid, ref, pdb_dir, out, timeout):
    outdir = Path(out) / pid
    t0 = time.time()
    cmd = [sys.executable, __file__, "--worker", str(Path(pdb_dir) / f"{pid}.pdb"), ref or "", str(outdir)]
    try:
        p = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)
    except subprocess.TimeoutExpired:
        return {"id": pid, "ref": ref, "status": "timeout", "seconds": timeout}
    sec = round(time.time() - t0, 1)
    if p.returncode == 0 and p.stdout.strip():
        d = json.loads(p.stdout.strip().splitlines()[-1])
        return {"id": pid, "ref": ref, "status": "ok", "seconds": sec, **d}
    tail = (p.stderr or p.stdout).strip().splitlines()[-6:]
    return {"id": pid, "ref": ref, "status": "failed", "seconds": sec, "error": "\n".join(tail)}


def main():
    if len(sys.argv) > 1 and sys.argv[1] == "--worker":
        worker(sys.argv[2], sys.argv[3] or None, sys.argv[4])
        return
    ap = argparse.ArgumentParser()
    ap.add_argument("--pdb-dir", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--jobs", type=int, default=4)
    ap.add_argument("--timeout", type=int, default=900)
    ap.add_argument("--only", nargs="*")
    ap.add_argument("--gate", action="store_true",
                    help="exit non-zero if any entry fails or loses protein residues")
    a = ap.parse_args()
    ids = a.only or list(SET)
    Path(a.out).mkdir(parents=True, exist_ok=True)
    with ThreadPoolExecutor(a.jobs) as ex:
        rows = list(ex.map(lambda i: run_one(i, SET[i], a.pdb_dir, a.out, a.timeout), ids))
    Path(a.out, "summary.json").write_text(json.dumps(rows, indent=2))
    md = ["| ID | ref | status | s | chains | atoms before→after | PDBQT atoms | warnings |", "|---|---|---|---|---|---|---|---|"]
    for r in rows:
        if r["status"] == "ok":
            md.append(f"| {r['id']} | {r['ref'] or '-'} | ok | {r['seconds']} | {','.join(r['chains'])} | "
                      f"{r['atoms_before']}→{r['atoms_after']} | {r['pdbqt_atoms']} | {r['n_warnings']} |")
        else:
            md.append(f"| {r['id']} | {r['ref'] or '-'} | **{r['status']}** | {r['seconds']} | | | | |")
    Path(a.out, "summary.md").write_text("\n".join(md) + "\n")
    print("\n".join(md))
    print(f"\n{sum(r['status']=='ok' for r in rows)}/{len(rows)} ok")
    problems = []
    for r in rows:
        if r["status"] != "ok":
            problems.append(f"{r['id']}: {r['status']}")
        else:
            if r["lost_residues_pdb"]:
                problems.append(f"{r['id']}: residues lost from prepared PDB: {r['lost_residues_pdb'][:5]}")
            if r["lost_residues_pdbqt"] and r["id"] not in EXPECTED_DEGRADED:
                problems.append(f"{r['id']}: residues lost from PDBQT: {r['lost_residues_pdbqt'][:5]}")
            if r["lost_residues_pdbqt"] and r["status_detail"] != "degraded":
                problems.append(f"{r['id']}: residues lost but the report does not say 'degraded'")
            if r["id"] in EXPECTED_DEGRADED and r["status_detail"] == "degraded":
                print(f"NOTE {r['id']} degraded as expected: {EXPECTED_DEGRADED[r['id']]} "
                      f"({r['n_dropped']} residues)")
    for p in problems:
        print("PROBLEM", p)
    if a.gate and problems:
        sys.exit(1)


if __name__ == "__main__":
    main()
