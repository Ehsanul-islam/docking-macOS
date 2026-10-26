"""Batch docking: run the Vina engine over many prepared ligands against one receptor and box.

Core budget (resources.cores()) is split sensibly: many ligands run one core each, side by side; a single ligand gets all
the cores. A ligand that fails or times out never stops the batch. Cancel kills the running engines and keeps what finished.
A run folder holds results.csv / results.json (ranked, best predicted binding affinity first), poses/, logs/, results/ (one
JSON per ligand, which makes a re-run resume where it stopped) and run.json (settings, engine version, cores).

Scores are Vina scores: a predicted binding affinity in kcal/mol for ranking, not a measured free energy.
"""
from __future__ import annotations

import csv
import hashlib
import json
import os
import re
import shutil
import signal
import subprocess
import time
from pathlib import Path

from . import __version__, cancel, resources
from .errors import Cancelled, InvalidInput, ToolFailure, UnsupportedFeature

DEFAULT_TIMEOUT_S = 3600
RESULT_RE = re.compile(r"REMARK VINA RESULT:\s+(-?\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)")


def find_engine(explicit: str | None = None) -> str:
    """Explicit path, else $DOCKCHEM_VINA, else a bundled engine/bin/vina, else `vina` on PATH."""
    candidates = [explicit, os.environ.get("DOCKCHEM_VINA"),
                  str(Path(__file__).resolve().parents[2] / "build" / "engine" / "vina"), shutil.which("vina")]
    for c in candidates:
        if c and Path(c).is_file() and os.access(c, os.X_OK):
            return str(c)
    raise UnsupportedFeature("the docking engine (AutoDock Vina) was not found. Install or bundle it, or set DOCKCHEM_VINA "
                             "to its path")


def engine_version(engine: str) -> str:
    try:
        out = subprocess.run([engine, "--version"], capture_output=True, text=True, timeout=20)
        return (out.stdout + out.stderr).strip().splitlines()[0]
    except Exception as e:  # noqa: BLE001
        raise ToolFailure(f"the docking engine did not run: {e}") from e


def plan_cores(n_ligands: int, cores: int) -> tuple[int, int]:
    """(simultaneous engine processes, threads each). Many ligands: one core each. One ligand: all cores."""
    n_ligands, cores = max(1, n_ligands), max(1, cores)
    procs = min(cores, n_ligands)
    return procs, max(1, cores // procs)


def parse_config(path) -> dict:
    """Read a vina.conf (key = value lines) such as quickprep writes."""
    out = {}
    for line in Path(path).read_text().splitlines():
        line = line.split("#")[0].strip()
        if "=" in line:
            k, v = (x.strip() for x in line.split("=", 1))
            out[k] = v
    return out


def _check_ligand(path: Path) -> None:
    text = path.read_text(errors="replace")
    if "ROOT" not in text or "ENDROOT" not in text or "TORSDOF" not in text:
        raise InvalidInput(f"{path.name} is not a ligand PDBQT (no ROOT/TORSDOF)")


def _safe(name: str) -> str:
    return "".join(c if c.isalnum() or c in "-_." else "_" for c in name) or "ligand"


def _settings(receptor: Path, center, size, exhaustiveness, num_modes, seed, scoring, engine_ver) -> dict:
    return {"receptor_sha1": hashlib.sha1(receptor.read_bytes()).hexdigest(), "center": [float(x) for x in center],
            "size": [float(x) for x in size], "exhaustiveness": int(exhaustiveness), "num_modes": int(num_modes),
            "seed": int(seed), "scoring": scoring, "engine": engine_ver}


def run_batch(receptor, ligands, outdir, center, size, exhaustiveness: int = 32, num_modes: int = 9, seed: int = 42,
              scoring: str = "vina", engine: str | None = None, timeout_s: float = DEFAULT_TIMEOUT_S,
              resume: bool = True, write_sdf: bool = True, progress=None) -> dict:
    receptor, outdir = Path(receptor), Path(outdir)
    if not receptor.is_file():
        raise InvalidInput(f"receptor PDBQT not found: {receptor}")
    ligs = [Path(l) for l in ligands]
    if not ligs:
        raise InvalidInput("no ligands to dock")
    if len(center) != 3 or len(size) != 3 or any(float(s) <= 0 for s in size):
        raise InvalidInput("center and size must each have 3 numbers, and every size must be positive")
    if not 1 <= int(exhaustiveness) <= 512:
        raise InvalidInput("exhaustiveness must be between 1 and 512")
    if not 1 <= int(num_modes) <= 50:
        raise InvalidInput("num_modes must be between 1 and 50")
    if scoring not in ("vina", "vinardo"):
        raise InvalidInput("scoring must be vina or vinardo (AutoDock4 scoring is not supported yet)")
    names, seen = [], {}
    for l in ligs:
        if not l.is_file():
            raise InvalidInput(f"ligand file not found: {l}")
        _check_ligand(l)
        n = _safe(l.stem)
        seen[n] = seen.get(n, 0) + 1
        names.append(n if seen[n] == 1 else f"{n}_{seen[n]}")
    eng = find_engine(engine)
    ever = engine_version(eng)
    cancel.check()

    for d in ("poses", "logs", "results"):
        (outdir / d).mkdir(parents=True, exist_ok=True)
    settings = _settings(receptor, center, size, exhaustiveness, num_modes, seed, scoring, ever)
    run_file = outdir / "run.json"
    done: dict[str, dict] = {}
    if run_file.exists():
        old = json.loads(run_file.read_text())
        if old.get("settings") != settings:
            if resume:
                raise InvalidInput(f"{outdir} holds a run with different settings (receptor, box, search or engine); "
                                   "use a new folder or set resume=False to overwrite it")
            for f in (outdir / "results").glob("*.json"):
                f.unlink()
        elif resume:
            for n in names:
                f = outdir / "results" / f"{n}.json"
                if f.exists():
                    done[n] = json.loads(f.read_text())
    cores = resources.cores()
    procs, threads = plan_cores(len([n for n in names if n not in done]), cores)
    run_file.write_text(json.dumps({"settings": settings, "cores": cores, "processes": procs, "threads_each": threads,
                                    "ligands": len(names), "started": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
                                    "app_version": __version__}, indent=1))

    todo = [(n, l) for n, l in zip(names, ligs) if n not in done]
    results = dict(done)
    running: dict = {}
    t_start, finished = time.time(), len(done)

    def launch(n, l):
        log = open(outdir / "logs" / f"{n}.log", "w")
        cmd = [eng, "--receptor", str(receptor), "--ligand", str(l), "--center_x", str(center[0]), "--center_y", str(center[1]),
               "--center_z", str(center[2]), "--size_x", str(size[0]), "--size_y", str(size[1]), "--size_z", str(size[2]),
               "--exhaustiveness", str(int(exhaustiveness)), "--num_modes", str(int(num_modes)), "--seed", str(int(seed)),
               "--cpu", str(threads), "--scoring", scoring, "--out", str(outdir / "poses" / f"{n}_out.pdbqt")]
        proc = subprocess.Popen(cmd, stdout=log, stderr=subprocess.STDOUT, start_new_session=True)
        running[n] = {"proc": proc, "log": log, "start": time.time(), "ligand": l}

    def kill(info):
        try:
            os.killpg(info["proc"].pid, signal.SIGKILL)
        except (ProcessLookupError, PermissionError):
            info["proc"].kill()
        info["proc"].wait()
        info["log"].close()

    def record(n, info, status, error=""):
        secs = round(time.time() - info["start"], 1)
        row = {"name": n, "status": status, "error": error, "seconds": secs, "affinity": None, "n_poses": 0,
               "poses": "", "scores": []}
        out = outdir / "poses" / f"{n}_out.pdbqt"
        if status == "docked":
            sc = [(float(a), float(b), float(c)) for a, b, c in RESULT_RE.findall(out.read_text())] if out.exists() else []
            if not sc:
                row.update(status="failed", error="the engine finished but wrote no poses")
            else:
                row.update(affinity=sc[0][0], n_poses=len(sc), poses=str(out), scores=[s[0] for s in sc])
                if write_sdf:
                    try:
                        from .formats import poses_to_sdf
                        sdf = outdir / "poses" / f"{n}_out.sdf"
                        poses_to_sdf(out, sdf)
                        row["sdf"] = str(sdf)
                    except Exception as e:  # noqa: BLE001
                        row["sdf_error"] = str(e)[:150]
        (outdir / "results" / f"{n}.json").write_text(json.dumps(row))
        results[n] = row

    queue = list(todo)
    try:
        while queue or running:
            tok = cancel.current()
            if tok is not None and tok.cancelled:
                for n, info in list(running.items()):
                    kill(info)
                _write_outputs(outdir, names, results, cores)
                raise Cancelled(f"job cancelled after {finished} of {len(names)} ligands; finished results are in {outdir}")
            while queue and len(running) < procs:
                launch(*queue.pop(0))
            time.sleep(0.15)
            for n, info in list(running.items()):
                rc = info["proc"].poll()
                if rc is None and time.time() - info["start"] > timeout_s:
                    kill(info)
                    record(n, info, "failed", f"exceeded the {timeout_s:.0f} s limit and was stopped")
                elif rc is not None:
                    info["log"].close()
                    if rc == 0:
                        record(n, info, "docked")
                    else:
                        tail = Path(info["log"].name).read_text(errors="replace").strip().splitlines()[-3:]
                        record(n, info, "failed", f"engine exited with code {rc}: " + " | ".join(tail)[:300])
                else:
                    continue
                running.pop(n)
                finished += 1
                if progress:
                    r = results[n]
                    eta = (time.time() - t_start) / max(1, finished - len(done)) * (len(names) - finished)
                    progress(finished / len(names), f"{n}: " + (f"{r['affinity']:.1f} kcal/mol" if r["affinity"] is not None
                                                                else "failed") + f" ({finished}/{len(names)}, ~{eta / 60:.0f} min left)")
    finally:
        for info in running.values():
            if info["proc"].poll() is None:
                kill(info)
    summary = _write_outputs(outdir, names, results, cores)
    summary["seconds"] = round(time.time() - t_start, 1)
    return summary


def _write_outputs(outdir: Path, names, results: dict, cores: int) -> dict:
    rows = [results[n] for n in names if n in results]
    docked = sorted((r for r in rows if r["status"] == "docked"), key=lambda r: (r["affinity"], r["name"]))
    ranked = [{**r, "rank": i + 1} for i, r in enumerate(docked)] + [r for r in rows if r["status"] != "docked"]
    with open(outdir / "results.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["rank", "name", "affinity_kcal_mol", "n_poses", "status", "seconds", "poses", "error"])
        for r in ranked:
            w.writerow([r.get("rank", ""), r["name"], r["affinity"] if r["affinity"] is not None else "", r["n_poses"],
                        r["status"], r["seconds"], r["poses"], r["error"]])
    summary = {"ligands": len(names), "docked": len(docked), "failed": sum(r["status"] == "failed" for r in rows),
               "not_run": len(names) - len(rows), "cores": cores, "results": ranked,
               "csv": str(outdir / "results.csv"), "outdir": str(outdir),
               "scores_are": "predicted binding affinity (Vina score, kcal/mol); lower is better; not a measured free energy"}
    (outdir / "results.json").write_text(json.dumps(summary, indent=1))
    return summary
