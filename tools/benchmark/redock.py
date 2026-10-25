#!/usr/bin/env python3
"""Redocking benchmark: how often does the full pipeline + unmodified Vina 1.2.7 put the known ligand back?

Per entry:
  receptor : prepared from the structure downloaded from the PDB, with the entry's ligand as site reference
  ligand   : the benchmark's generated start conformer (NOT the crystal pose), prepared by our ligand pipeline
  box      : centred on the crystal ligand, ligand extent + 2 x PAD on each axis (min MIN_BOX)
  docking  : vina 1.2.7, fixed seed, num_modes 9
  score    : heavy-atom symmetry-aware RMSD of each pose to the crystal ligand; success = top pose <= 2.0 A

  python tools/benchmark/redock.py --zip pb.zip --data /tmp/pdb300 --vina /path/to/vina --out /tmp/redock --set astex
  python tools/benchmark/redock.py --out /tmp/redock --report
"""
import argparse
import collections
import json
import re
import shutil
import subprocess
import sys
import tempfile
import time
import zipfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "batch"))
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "chemhelper"))

PAD = 6.0          # A of box around the crystal ligand on every side
MIN_BOX = 20.0
SUCCESS = 2.0


def read_sdf_from_zip(zpath, name):
    with zipfile.ZipFile(zpath) as z:
        return z.read(name).decode()


def nearest_copy(path, code, ref) -> str:
    """chain:resname:seq of the copy of the ligand in the structure that sits at the benchmark's crystal pose.
    Entries with several copies (chains A and B) would otherwise be built around the wrong one."""
    import math

    import gemmi
    c = ref.GetConformer().GetPositions().mean(axis=0)
    st = gemmi.read_structure(str(path))
    best = (1e9, None)
    for ch in st[0]:
        for r in ch:
            if r.name == code and r.het_flag == "H":
                p = [a.pos for a in r if a.element.name != "H"]
                if p:
                    cen = (sum(x.x for x in p) / len(p), sum(x.y for x in p) / len(p), sum(x.z for x in p) / len(p))
                    d = math.dist(cen, c)
                    if d < best[0]:
                        best = (d, f"{ch.name}:{r.name}:{r.seqid.num}")
    return best[1] or code


def strip_h(m):
    """Remove every hydrogen. RDKit's RemoveHs keeps ones that define double-bond stereo ('[H]/N=C'), which the
    benchmark's start conformers contain."""
    from rdkit import Chem
    rw = Chem.RWMol(m)
    for i in sorted((a.GetIdx() for a in rw.GetAtoms() if a.GetAtomicNum() == 1), reverse=True):
        rw.RemoveAtom(i)
    out = rw.GetMol()
    out.UpdatePropertyCache(strict=False)
    return out


def heavy_rmsd(pose, ref):
    """Symmetry-aware heavy-atom RMSD without alignment. Falls back to a bond-order-blind match when the two
    molecules differ in charge or aromaticity assignment (same skeleton)."""
    from rdkit import Chem
    from rdkit.Chem import rdMolAlign
    try:
        return rdMolAlign.CalcRMS(strip_h(pose), strip_h(ref))
    except Exception:  # noqa: BLE001
        pass

    def skeleton(m):
        m = Chem.RWMol(strip_h(m))
        for b in m.GetBonds():
            b.SetBondType(Chem.BondType.SINGLE)
            b.SetIsAromatic(False)
        for a in m.GetAtoms():
            a.SetFormalCharge(0)
            a.SetIsAromatic(False)
            a.SetNoImplicit(True)
            a.SetNumExplicitHs(0)
        return m.GetMol()
    sp, sr = skeleton(pose), skeleton(ref)
    matches = sr.GetSubstructMatches(sp, uniquify=False, maxMatches=20000, useChirality=False)
    if not matches:
        raise ValueError("no skeleton match")
    import numpy as np
    P = sp.GetConformer().GetPositions()
    R = sr.GetConformer().GetPositions()
    return min(float(np.sqrt(((P - R[list(m)]) ** 2).sum(axis=1).mean())) for m in matches)


def work(item):
    from rdkit import Chem, RDLogger
    from rdkit.Chem import rdMolAlign
    RDLogger.DisableLog("rdApp.*")
    from dockchem.errors import HelperError
    from dockchem.formats import read_pdbqt_poses
    from dockchem.ligprep import LigandPrepOptions, prepare_ligand
    from dockchem.recprep import ReceptorPrepOptions, prepare_receptor

    key = f"{item['id']}_{item['ligand']}"
    base = f"{item['set']}/{key}/{key}"
    out = {"stage": "start", "set": item["set"]}
    ref = Chem.MolFromMolBlock(read_sdf_from_zip(item["zip"], base + "_ligand.sdf"), removeHs=False)
    start = Chem.MolFromMolBlock(read_sdf_from_zip(item["zip"], base + "_ligand_start_conf.sdf"), removeHs=False)
    if ref is None or start is None:
        return {**out, "category": "benchmark_file_unreadable", "error": "SDF did not parse"}
    work = Path(tempfile.mkdtemp(prefix="redock_"))
    try:
        out["stage"] = "receptor"
        t0 = time.time()
        if item.get("baseline"):
            # Standard route on the benchmark's own cleaned protein: Meeko straight on the file, nothing else.
            (work / "rec").mkdir()
            prot = work / "rec" / "protein.pdb"
            prot.write_text(read_sdf_from_zip(item["zip"], base + "_protein.pdb"))
            env = {**__import__("os").environ, "http_proxy": "http://127.0.0.1:9", "https_proxy": "http://127.0.0.1:9"}
            exe = str(Path(sys.executable).parent / "mk_prepare_receptor.py")
            exe = exe if Path(exe).exists() else "mk_prepare_receptor.py"
            m = subprocess.run([exe, "-i", str(prot), "-o", str(work / "rec" / "receptor"), "-p", "-a"],
                               capture_output=True, text=True, env=env, timeout=900)
            pq = work / "rec" / "receptor.pdbqt"
            if m.returncode != 0 or not pq.exists():
                return {**out, "category": "receptor_failed", "error": "Meeko: " + (m.stdout + m.stderr)[-200:].replace("\n", " ")}
            from types import SimpleNamespace
            rec = SimpleNamespace(pdbqt=pq, report={"status": "baseline", "dropped_residues": []})
        else:
            path = next(Path(item["data"]).glob(item["id"] + ".*"))
            spec = nearest_copy(path, item["ligand"], ref)          # the copy the benchmark's coordinates belong to
            out["reference"] = spec
            try:
                rec = prepare_receptor(path, work / "rec", ReceptorPrepOptions(reference_ligand=spec))
            except HelperError as e:
                return {**out, "category": "receptor_failed", "error": f"{e.kind}: {e.message[:200]}"}
        out["receptor_status"] = rec.report["status"]
        out["n_dropped"] = len(rec.report["dropped_residues"])
        out["receptor_seconds"] = round(time.time() - t0, 1)

        out["stage"] = "ligand"
        lig = [p for p in prepare_ligand(start, key, LigandPrepOptions(protonate=False)) if p.status != "failed"]
        if not lig:
            return {**out, "category": "ligand_failed", "error": "ligand preparation failed"}
        (work / "lig.pdbqt").write_text(lig[0].pdbqt)

        out["stage"] = "docking"
        xyz = ref.GetConformer().GetPositions()
        heavy = [i for i, a in enumerate(ref.GetAtoms()) if a.GetAtomicNum() > 1]
        c = xyz[heavy].mean(axis=0)
        ext = xyz[heavy].max(axis=0) - xyz[heavy].min(axis=0)
        size = [max(MIN_BOX, float(e) + 2 * PAD) for e in ext]
        cmd = [item["vina"], "--receptor", str(rec.pdbqt), "--ligand", str(work / "lig.pdbqt"),
               "--center_x", f"{c[0]:.3f}", "--center_y", f"{c[1]:.3f}", "--center_z", f"{c[2]:.3f}",
               "--size_x", f"{size[0]:.1f}", "--size_y", f"{size[1]:.1f}", "--size_z", f"{size[2]:.1f}",
               "--seed", "42", "--num_modes", "9", "--cpu", str(item.get("cpu", 1)),
               "--out", str(work / "out.pdbqt")]
        t0 = time.time()
        r = subprocess.run(cmd, capture_output=True, text=True, timeout=item.get("vina_timeout", 1500))
        out["vina_seconds"] = round(time.time() - t0, 1)
        if r.returncode != 0 or not (work / "out.pdbqt").exists():
            return {**out, "category": "vina_failed", "error": (r.stdout + r.stderr)[-300:].replace("\n", " ")}
        scores = [float(m.group(1)) for m in re.finditer(r"^\s+\d+\s+(-?\d+\.\d+)", r.stdout, re.M)]

        out["stage"] = "scoring"
        poses, _ = read_pdbqt_poses(work / "out.pdbqt")
        rmsds = []
        for pose in poses:
            try:
                rmsds.append(round(heavy_rmsd(pose, ref), 2))
            except Exception:  # noqa: BLE001
                rmsds.append(None)
        good = [x for x in rmsds if x is not None]
        if not good:
            return {**out, "category": "rmsd_failed", "error": "pose and crystal ligand could not be matched"}
        out.update(category="success" if rmsds[0] is not None and rmsds[0] <= SUCCESS else "miss",
                   top1_rmsd=rmsds[0], best_rmsd=min(good), rmsds=rmsds, scores=scores[:len(rmsds)],
                   best_in_top3=min(x for x in rmsds[:3] if x is not None) if any(x is not None for x in rmsds[:3]) else None,
                   n_torsions=lig[0].pdbqt.count("BRANCH") // 2, heavy_atoms=len(heavy),
                   ligand_warnings=lig[0].warnings[:3])
        return out
    except subprocess.TimeoutExpired:
        return {**out, "category": "vina_timeout", "error": "vina exceeded its time limit"}
    finally:
        shutil.rmtree(work, ignore_errors=True)


def report(path):
    rows = [json.loads(l) for l in Path(path).read_text().splitlines() if l.strip()]
    n = len(rows)
    cat = collections.Counter(r["category"] for r in rows)
    print(f"\n{n} entries")
    for c, k in cat.most_common():
        print(f"  {c:26s}{k:5d}  {100 * k / n:5.1f}%")
    scored = [r for r in rows if r["category"] in ("success", "miss")]
    if scored:
        s1 = sum(r["top1_rmsd"] <= SUCCESS for r in scored)
        s3 = sum((r["best_in_top3"] or 99) <= SUCCESS for r in scored)
        s9 = sum(r["best_rmsd"] <= SUCCESS for r in scored)
        print(f"\ndocked and scored: {len(scored)}")
        print(f"  top-1 within {SUCCESS} A : {s1}/{len(scored)} = {100 * s1 / len(scored):.1f}%   (of all {n} attempted: {100 * s1 / n:.1f}%)")
        print(f"  best of top 3          : {s3}/{len(scored)} = {100 * s3 / len(scored):.1f}%")
        print(f"  best of all 9 poses    : {s9}/{len(scored)} = {100 * s9 / len(scored):.1f}%  (search succeeded even if the scoring ranked it lower)")
        deg = [r for r in scored if r.get("receptor_status") == "degraded"]
        if deg:
            print(f"  entries with a degraded receptor: {len(deg)}, top-1 success {sum(r['top1_rmsd'] <= SUCCESS for r in deg)}/{len(deg)}")
        bins = collections.defaultdict(lambda: [0, 0])
        for r in scored:
            b = "0-5 torsions" if r["n_torsions"] <= 5 else "6-10 torsions" if r["n_torsions"] <= 10 else ">10 torsions"
            bins[b][0] += 1
            bins[b][1] += r["top1_rmsd"] <= SUCCESS
        for b, (t, s) in sorted(bins.items()):
            print(f"  {b:14s}: {s}/{t} = {100 * s / t:.0f}%")
    errs = collections.defaultdict(list)
    for r in rows:
        if r["category"] not in ("success", "miss"):
            errs[(r["category"], re.sub(r"\d+", "N", r.get("error", ""))[:90])].append(r["id"])
    for (c, e), ids in sorted(errs.items(), key=lambda kv: -len(kv[1])):
        print(f"  {len(ids):3d} {c} / {e}   e.g. {', '.join(ids[:4])}")
    t = sorted(r.get("vina_seconds", 0) for r in scored)
    if t:
        print(f"vina seconds per entry: median {t[len(t) // 2]:.0f}, max {t[-1]:.0f}")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--zip")
    ap.add_argument("--data")
    ap.add_argument("--vina")
    ap.add_argument("--out", required=True)
    ap.add_argument("--set", default="astex", choices=["astex", "posebusters", "both"])
    ap.add_argument("--manifest", default=str(Path(__file__).resolve().parents[1] / "batch" / "sets" / "widely_used_300.json"))
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--timeout", type=float, default=3000)
    ap.add_argument("--baseline", action="store_true", help="standard route: benchmark protein file straight through Meeko")
    ap.add_argument("--report", action="store_true")
    a = ap.parse_args()
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    if not a.report:
        from runner import run_batch
        items = json.load(open(a.manifest))["items"]
        names = {"astex": ["astex_diverse_set"], "posebusters": ["posebusters_benchmark_set"],
                 "both": ["astex_diverse_set", "posebusters_benchmark_set"]}[a.set]
        items = [dict(i, set=i["source"], zip=a.zip, data=a.data, vina=a.vina, cpu=1, baseline=a.baseline) for i in items if i["source"] in names]
        if a.limit:
            items = items[:a.limit]
        print(run_batch(items, "redock:work", out / "redock.jsonl", a.workers, a.timeout, progress_every=5))
    if (out / "redock.jsonl").exists():
        report(out / "redock.jsonl")
