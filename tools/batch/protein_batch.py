#!/usr/bin/env python3
"""Run a set of PDB entries through assessment + receptor preparation and summarise what happens.

  python tools/batch/protein_batch.py --manifest tools/batch/sets/widely_used_300.json --data /tmp/pdb --download
  python tools/batch/protein_batch.py --manifest ... --data /tmp/pdb --out /tmp/prot_run --workers 3
  python tools/batch/protein_batch.py --out /tmp/prot_run --report

Per entry: the assessment verdict (ready / fixable / unusable), then, unless unusable, the full preparation with
the entry's ligand as the site reference. Outcomes: ok, degraded (residues left out, named), unusable (refused with
reasons), failed (an error, grouped by message), timeout, crash. The report also counts how often the built-in
PDBQT writer and Meeko agree.
"""
import argparse
import collections
import json
import re
import shutil
import sys
import tempfile
import time
import urllib.request
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "chemhelper"))

DATA_DIR = None


def download(items, data: Path, log=print):
    data.mkdir(parents=True, exist_ok=True)
    for it in items:
        pid = it["id"]
        if any((data / f"{pid}{ext}").exists() and (data / f"{pid}{ext}").stat().st_size > 1000 for ext in (".pdb", ".cif")):
            continue
        for ext in (".pdb", ".cif"):          # large structures have no legacy PDB file
            url = f"https://files.rcsb.org/download/{pid}{ext}"
            for attempt in range(3):
                try:
                    with urllib.request.urlopen(url, timeout=60) as r:
                        body = r.read()
                    if len(body) > 1000:
                        (data / f"{pid}{ext}").write_bytes(body)
                        break
                except Exception:  # noqa: BLE001
                    time.sleep(1 + attempt)
            if (data / f"{pid}{ext}").exists():
                break
        else:
            log(f"  could not download {pid}")
    n = sum(1 for it in items if any((data / f"{it['id']}{e}").exists() for e in (".pdb", ".cif")))
    log(f"{n}/{len(items)} structures available in {data}")


def work(item):
    """Assess and prepare one entry (runs in a worker process)."""
    from rdkit import RDLogger
    RDLogger.DisableLog("rdApp.*")
    from dockchem.assess import assess
    from dockchem.errors import HelperError, UnusableStructure
    from dockchem.recprep import ReceptorPrepOptions, prepare_receptor
    path = next((Path(item["data"]) / f"{item['id']}{e}" for e in (".pdb", ".cif")
                 if (Path(item["data"]) / f"{item['id']}{e}").exists()), None)
    if path is None:
        return {"category": "missing_file", "error": "structure not downloaded"}
    out = {"stats": {}, **{k: item.get(k) for k in ("organism_group", "method", "uniprot", "structures_for_protein",
                                                     "n_atoms", "ligand")}}
    a = assess(path)
    out["verdict"] = a["verdict"]
    out["finding_codes"] = sorted({f["code"] for f in a["findings"]})
    out["stats"] = {k: a["stats"].get(k) for k in ("atoms", "protein_residues", "distorted_residues", "chain_breaks",
                                                    "residues_with_missing_atoms", "method", "resolution")}
    if a["verdict"] == "unusable":
        out.update(category="unusable", error=a["summary"][:200])
        return out
    work_dir = Path(tempfile.mkdtemp(prefix="prot_"))
    try:
        t0 = time.time()
        try:
            res = prepare_receptor(path, work_dir, ReceptorPrepOptions(reference_ligand=item.get("ligand")))
        except UnusableStructure as e:
            out.update(category="unusable", error=str(e)[:200])
            return out
        except HelperError as e:
            out.update(category="failed", error_kind=e.kind, error=e.message[:300], log_tail=(e.log or "")[-400:])
            return out
        rp = res.report
        cc = (rp.get("pdbqt") or {}).get("crosscheck") or {}
        out.update(category="ok" if rp["status"] == "ok" else "degraded",
                   dropped=[d["residue"] + " :: " + d["reason"][:60] for d in rp["dropped_residues"]][:20],
                   n_dropped=len(rp["dropped_residues"]), n_warnings=len(rp["warnings"]),
                   writer=(rp.get("pdbqt") or {}).get("writer"), crosscheck=cc.get("result"),
                   crosscheck_examples=cc.get("examples", [])[:3] if cc.get("result") == "disagree" else [],
                   crosscheck_reason=cc.get("reason", "")[:150] if cc.get("result") == "unavailable" else "",
                   n_flips=sum(1 for c in rp["changes"] if c["stage"] == "flip"),
                   n_ligands=len(res.ligands),
                   ligand_bond_orders=dict(collections.Counter((l.get("bond_orders") or {}).get("status", "?") for l in res.ligands)),
                   warnings=rp["warnings"][:25], assessment_findings=[f["code"] for f in a["findings"] if f["severity"] in ("fix", "block")],
                   prep_seconds=round(time.time() - t0, 1))
        return out
    finally:
        shutil.rmtree(work_dir, ignore_errors=True)


def norm(msg):
    return re.sub(r"\d+", "N", (msg or "").replace("\n", " "))[:90]


def report(path):
    rows = [json.loads(l) for l in Path(path).read_text().splitlines() if l.strip()]
    n = len(rows)
    print(f"\n{n} entries")
    cat = collections.Counter(r["category"] for r in rows)
    for c, k in cat.most_common():
        print(f"  {c:12s}{k:6d}  {100 * k / n:5.1f}%")
    print("\nassessment verdicts:", dict(collections.Counter(r.get("verdict", "?") for r in rows)))
    prepared = [r for r in rows if r["category"] in ("ok", "degraded")]
    print(f"prepared (ok + degraded): {len(prepared)}/{n}")
    cc = collections.Counter(r.get("crosscheck") for r in prepared)
    print("PDBQT writer:", dict(collections.Counter(r.get("writer") for r in prepared)), "| Meeko cross-check:", dict(cc))
    t = sorted(r.get("seconds", 0) for r in rows)
    if t:
        print(f"seconds per entry: median {t[len(t) // 2]:.0f}, p95 {t[int(len(t) * .95)]:.0f}, max {t[-1]:.0f}")
    groups = collections.defaultdict(list)
    for r in rows:
        if r["category"] in ("failed", "timeout", "crash", "exception", "missing_file"):
            groups[(r["category"], norm(r.get("error")))].append(r["id"])
    if groups:
        print("\nfailure groups:")
        for (c, e), ids in sorted(groups.items(), key=lambda kv: -len(kv[1])):
            print(f"  {len(ids):4d}  {c} / {e}   e.g. {', '.join(ids[:5])}")
    unus = collections.defaultdict(list)
    for r in rows:
        if r["category"] == "unusable":
            unus[norm(r.get("error"))[:70]].append(r["id"])
    if unus:
        print("\nunusable (refused with reasons):")
        for e, ids in sorted(unus.items(), key=lambda kv: -len(kv[1])):
            print(f"  {len(ids):4d}  {e}   e.g. {', '.join(ids[:5])}")
    deg = collections.Counter()
    for r in prepared:
        for d in r.get("dropped", []):
            deg[d.split(" :: ")[1][:60]] += 1
    if deg:
        print("\nresidues left out, by reason:")
        for e, k in deg.most_common(10):
            print(f"  {k:5d}  {e}")
    dis = [r for r in prepared if r.get("crosscheck") == "disagree"]
    if dis:
        print("\nwriter/Meeko disagreements:")
        for r in dis[:8]:
            print("  ", r["id"], r["crosscheck_examples"][:1])


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--manifest")
    ap.add_argument("--data", help="directory holding the downloaded structures")
    ap.add_argument("--out")
    ap.add_argument("--download", action="store_true")
    ap.add_argument("--workers", type=int, default=3)
    ap.add_argument("--timeout", type=float, default=900)
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--report", action="store_true")
    a = ap.parse_args()
    out = Path(a.out) if a.out else None
    items = json.load(open(a.manifest))["items"] if a.manifest else []
    if a.limit:
        items = items[:a.limit]
    if a.download:
        download(items, Path(a.data))
    if out and not a.report and items:
        out.mkdir(parents=True, exist_ok=True)
        from runner import run_batch
        for it in items:
            it["data"] = str(a.data)
        print(run_batch(items, "protein_batch:work", out / "proteins.jsonl", a.workers, a.timeout, progress_every=10))
    if out and (out / "proteins.jsonl").exists():
        report(out / "proteins.jsonl")
