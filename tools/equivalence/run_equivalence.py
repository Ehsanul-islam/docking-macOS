#!/usr/bin/env python3
"""Engine-equivalence harness (plan §12.1).

Runs a reference and a candidate Vina binary on the same cases with fixed seeds
and compares the output poses and scores.

  same architecture:   output files must be byte-identical (--mode exact)
  cross architecture:  top score within --tol kcal/mol        (--mode tolerance)
"""
import argparse
import json
import re
import subprocess
import sys
import tempfile
from pathlib import Path

SCORE_RE = re.compile(r"^REMARK VINA RESULT:\s+(-?\d+\.\d+)", re.M)


def run_vina(binary, case, outdir, base):
    out = Path(outdir) / f"{case['name']}.pdbqt"
    cmd = [binary, "--receptor", str(base / case["receptor"]), "--ligand", str(base / case["ligand"]),
           "--center_x", str(case["center"][0]), "--center_y", str(case["center"][1]),
           "--center_z", str(case["center"][2]),
           "--size_x", str(case["size"][0]), "--size_y", str(case["size"][1]), "--size_z", str(case["size"][2]),
           "--scoring", case.get("scoring", "vina"), "--seed", str(case["seed"]),
           "--exhaustiveness", str(case.get("exhaustiveness", 8)), "--out", str(out)]
    subprocess.run(cmd, check=True, capture_output=True)
    return out.read_text()


def scores(text):
    return [float(s) for s in SCORE_RE.findall(text)]


def compare(ref_text, cand_text, mode, tol):
    if mode == "exact":
        return ref_text == cand_text, "outputs differ" if ref_text != cand_text else ""
    r, c = scores(ref_text), scores(cand_text)
    if not r or not c:
        return False, "no scores parsed"
    d = abs(r[0] - c[0])
    return d <= tol, f"top-score delta {d:.3f} > {tol}" if d > tol else ""


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--reference", required=True)
    ap.add_argument("--candidate", required=True)
    ap.add_argument("--cases", required=True)
    ap.add_argument("--mode", choices=["exact", "tolerance"], default="exact")
    ap.add_argument("--tol", type=float, default=0.1)
    a = ap.parse_args(argv)
    cases_path = Path(a.cases).resolve()
    base = cases_path.parent
    failed = 0
    with tempfile.TemporaryDirectory() as td:
        for case in json.loads(cases_path.read_text())["cases"]:
            (Path(td) / "ref").mkdir(exist_ok=True); (Path(td) / "cand").mkdir(exist_ok=True)
            ref = run_vina(a.reference, case, Path(td) / "ref", base)
            cand = run_vina(a.candidate, case, Path(td) / "cand", base)
            ok, why = compare(ref, cand, a.mode, a.tol)
            print(f"{'PASS' if ok else 'FAIL'} {case['name']} {why}")
            failed += not ok
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
