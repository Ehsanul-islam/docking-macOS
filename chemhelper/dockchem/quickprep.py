"""One-panel quick prep (gap 6): a receptor and a ligand in, a docking-ready folder out.

  quick_prep(receptor, ligand, outdir, preset="standard", site=None)

Everything the separate steps would ask for is decided here with defaults a non-expert can trust, and every decision
is returned so the UI can show it: structure verdict, binding site (and why), box, receptor/ligand status, engine
settings, warnings. Nothing is hidden and nothing is guessed silently: a structure that cannot be docked into is
refused with the reasons, no site means blind docking and says so.

`ligand`: a SMILES string, or a path to SDF / MOL / SMI / PDB (several molecules = a library).
`site`  : None (largest drug-like ligand in the structure), {"ligand": "A:OHT:600"} or {"ligand": "OHT"},
          {"center": [x, y, z], "size": [sx, sy, sz]}, or {"residues": ["A:394", "A:353"]} (box around them).
"""
from __future__ import annotations

import json
import time
from pathlib import Path

import gemmi
import numpy as np

from . import hetatm
from .assess import assess
from .errors import InvalidInput, UnusableStructure
from .ligprep import LigandPrepOptions, prepare_library, prepare_ligand, read_records
from .recprep import ReceptorPrepOptions, prepare_receptor

PRESETS = {
    "fast": dict(exhaustiveness=8, num_modes=9, label="Fast (exhaustiveness 8): a quick look"),
    "standard": dict(exhaustiveness=32, num_modes=9, label="Standard (exhaustiveness 32): recommended"),
    "thorough": dict(exhaustiveness=64, num_modes=20, label="Thorough (exhaustiveness 64): flexible or large sites"),
}
PADDING = 6.0           # A around the reference ligand on every side
MIN_BOX = 20.0
BIG_BOX = 40.0          # above this per side, accuracy drops unless the search is lengthened


def _largest_ligand(path: Path):
    st = gemmi.read_structure(str(path))
    best, others = (0, None), []
    for ch in st[0]:
        for r in ch:
            if r.het_flag != "H":
                continue
            heavy = sum(1 for a in r if not a.element.is_hydrogen)
            if hetatm.classify(r.name, heavy) == "ligand":
                others.append(r.name)
                if heavy > best[0]:
                    best = (heavy, f"{ch.name}:{r.name}:{r.seqid.num}")
    return best[1], sorted({n for n in others if best[1] and n != best[1].split(":")[1]})


def _box_from_points(pts: np.ndarray, padding: float):
    lo, hi = pts.min(0), pts.max(0)
    center = (lo + hi) / 2
    size = np.maximum(hi - lo + 2 * padding, MIN_BOX)
    return [round(float(x), 3) for x in center], [round(float(x), 1) for x in size]


def _points_of_residues(path: Path, keys: list[str]) -> np.ndarray:
    st = gemmi.read_structure(str(path))
    want = set()
    for k in keys:
        chain, _, num = str(k).partition(":")
        if not num:
            raise InvalidInput(f"residue '{k}' must look like 'A:394'")
        want.add((chain, num.strip()))
    pts = [[a.pos.x, a.pos.y, a.pos.z] for ch in st[0] for r in ch
           if (ch.name, f"{r.seqid.num}{r.seqid.icode.strip()}") in want for a in r if not a.element.is_hydrogen]
    if not pts:
        raise InvalidInput(f"none of the residues {sorted(k for _, k in want)} were found in the structure")
    return np.array(pts)


def quick_prep(receptor, ligand, outdir, preset: str = "standard", site: dict | None = None, seed: int = 42,
               allow_unusable: bool = False, receptor_options: dict | None = None,
               ligand_options: dict | None = None) -> dict:
    if preset not in PRESETS:
        raise InvalidInput(f"preset must be one of {', '.join(PRESETS)}, not {preset!r}")
    t0 = time.time()
    receptor, outdir = Path(receptor), Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    out: dict = {"preset": preset, "steps": [], "warnings": []}

    def step(name, t):
        out["steps"].append({"step": name, "seconds": round(time.time() - t, 1)})

    # 1. is this structure usable?
    t = time.time()
    a = assess(receptor)
    out["assessment"] = {"verdict": a["verdict"], "summary": a["summary"]}
    if a["verdict"] == "unusable" and not allow_unusable:
        raise UnusableStructure("this structure cannot be docked into: " + a["summary"], a)
    step("assess", t)

    # 2. where to dock
    t = time.time()
    site = dict(site or {})
    ref_spec, box_pts, how = None, None, ""
    if "center" in site:
        c = [float(x) for x in site["center"]]
        s = [float(x) for x in site.get("size", [22, 22, 22])]
        center, size, how = c, s, "box given by the user"
    elif "residues" in site:
        center, size = _box_from_points(_points_of_residues(receptor, site["residues"]), PADDING)
        how = "box around residues " + ", ".join(map(str, site["residues"]))
    else:
        ref_spec = site.get("ligand")
        if ref_spec is None:
            ref_spec, others = _largest_ligand(receptor)
            how = f"around the largest ligand in the structure ({ref_spec})" if ref_spec else ""
            if others:
                out["warnings"].append("other ligands also present: " + ", ".join(others) + "; use site={'ligand': ...} to pick one")
        else:
            how = f"around the ligand {ref_spec}"
    step("site", t)

    # 3. receptor
    t = time.time()
    ropts = ReceptorPrepOptions(reference_ligand=ref_spec, allow_unusable=allow_unusable, **(receptor_options or {}))
    rec = prepare_receptor(receptor, outdir / "receptor", ropts)
    out["receptor"] = {"status": rec.report["status"], "pdbqt": str(rec.pdbqt), "pdb": str(rec.pdb),
                       "left_out": [d["residue"] + ": " + d["reason"] for d in rec.report["dropped_residues"]],
                       "report": str(outdir / "receptor" / "prep_report.json")}
    (outdir / "receptor" / "prep_report.json").write_text(json.dumps(rec.report, indent=1, default=str))
    out["warnings"] += rec.report["warnings"][:20]
    step("receptor", t)

    # box from the reference ligand's coordinates, or blind docking
    if ref_spec and rec.ligands and "center" not in site and "residues" not in site:
        pts = []
        for l in rec.ligands:
            pts += [[at.pos.x, at.pos.y, at.pos.z] for r in gemmi.read_structure(l["path"])[0][0] for at in r
                    if not at.element.is_hydrogen]
        center, size = _box_from_points(np.array(pts), PADDING)
    elif "center" not in site and "residues" not in site:
        st = gemmi.read_structure(str(rec.pdb))
        pts = np.array([[a.pos.x, a.pos.y, a.pos.z] for ch in st[0] for r in ch for a in r if not a.element.is_hydrogen])
        center, size = _box_from_points(pts, 0.0)
        how = "no binding site given or found: blind docking over the whole protein"
        out["warnings"].append("blind docking: accuracy is much lower than docking into a known site, and runs take longer")
    out["site"] = {"how": how, "center": center, "size": size}
    if max(size) > BIG_BOX:
        out["warnings"].append(f"the search box is {max(size):.0f} A on its longest side; use the Thorough preset or a smaller site")

    # 4. ligand(s)
    t = time.time()
    lopts = LigandPrepOptions(**(ligand_options or {}))
    ldir = outdir / "ligands"
    ldir.mkdir(exist_ok=True)
    if isinstance(ligand, (str, Path)) and Path(str(ligand)).exists() and Path(str(ligand)).is_file():
        lib = prepare_library(read_records(ligand), lopts)
    elif isinstance(ligand, str) and ligand.strip():
        lib = {"prepared": [p for p in prepare_ligand(ligand, "ligand", lopts)], "rows": []}
        lib["rows"] = [p.summary() for p in lib["prepared"]]
    else:
        raise InvalidInput("ligand must be a SMILES string or the path of an SDF / MOL / SMI / PDB file")
    files = []
    for p in lib["prepared"]:
        if p.status != "failed" and p.pdbqt:
            f = ldir / f"{''.join(c if c.isalnum() or c in '-_.' else '_' for c in p.name)}.pdbqt"
            f.write_text(p.pdbqt)
            files.append(str(f))
    failed = [r for r in lib["rows"] if r["status"] == "failed"]
    out["ligands"] = {"prepared": len(files), "failed": [{"name": r["name"], "error": r["error"]} for r in failed],
                      "files": files}
    if not files:
        raise InvalidInput("no ligand could be prepared: " + "; ".join(f"{r['name']}: {r['error']}" for r in failed[:3]))
    step("ligands", t)

    # 5. engine settings
    pr = PRESETS[preset]
    conf = outdir / "vina.conf"
    conf.write_text("\n".join([f"receptor = {rec.pdbqt}", f"center_x = {center[0]}", f"center_y = {center[1]}",
                               f"center_z = {center[2]}", f"size_x = {size[0]}", f"size_y = {size[1]}", f"size_z = {size[2]}",
                               f"exhaustiveness = {pr['exhaustiveness']}", f"num_modes = {pr['num_modes']}",
                               f"seed = {seed}", "cpu = 0"]) + "\n")
    out["engine"] = {"scoring": "vina", **{k: pr[k] for k in ("exhaustiveness", "num_modes")}, "seed": seed,
                     "config": str(conf), "preset_label": pr["label"]}
    out["ready_to_dock"] = True
    out["seconds"] = round(time.time() - t0, 1)
    (outdir / "quickprep_summary.json").write_text(json.dumps(out, indent=1, default=str))
    return out
