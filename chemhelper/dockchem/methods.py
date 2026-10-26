"""RPC method registry. Phase 0 exposes system + project + a ligand smoke method."""
import json
import platform
import sys

from . import __version__
from .errors import InvalidInput
from .project import Project
from .provenance import record


def ping(params, srv):
    return {"result": "pong"}


def versions(params, srv):
    import gemmi
    import rdkit
    return {"result": {
        "dockchem": __version__,
        "python": sys.version.split()[0],
        "rdkit": rdkit.__version__,
        "gemmi": gemmi.__version__,
        "arch": platform.machine(),
    }}


def _req(params, key):
    """A required parameter, or InvalidInput (never a bare KeyError)."""
    v = params.get(key)
    if v is None or (isinstance(v, str) and not v.strip()):
        raise InvalidInput(f"missing required parameter '{key}'")
    return v


def _text(params, key):
    v = _req(params, key)
    if not isinstance(v, str):
        raise InvalidInput(f"parameter '{key}' must be text, not {type(v).__name__}")
    try:
        v.encode("utf-8")
    except UnicodeEncodeError as e:
        raise InvalidInput(f"parameter '{key}' contains characters that are not valid text") from e
    return v


def _path(params, key="path"):
    """A required input file that exists, is readable and is not absurdly large."""
    from .limits import check_input_file, check_structure_text
    v = _req(params, key)
    if not isinstance(v, str):
        raise InvalidInput(f"parameter '{key}' must be a file path (text), not {type(v).__name__}")
    p = check_input_file(v, what=f"'{key}'")
    check_structure_text(p, what=f"'{key}'")
    return str(p)


def sleep_job(params, srv):
    """Cancellable test job exercising progress notifications."""
    import time
    from .errors import Cancelled
    job_id = params.get("job_id")
    steps = int(params.get("steps", 5))
    for i in range(steps):
        if job_id in srv.cancelled:
            raise Cancelled("job cancelled")
        srv.progress(job_id, (i + 1) / steps)
        time.sleep(float(params.get("interval", 0.01)))
    return {"result": "done"}


def ligand_inspect(params, srv):
    """Parse a SMILES and report basic properties (smoke test for RDKit)."""
    from rdkit import Chem
    from rdkit.Chem import Descriptors, rdMolDescriptors
    smi = params.get("smiles")
    smi = _text(params, "smiles") if smi else smi
    if not smi:
        raise InvalidInput("smiles required (a text string)")
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        raise InvalidInput(f"cannot parse SMILES: {smi}")
    return {
        "result": {
            "inchikey": Chem.MolToInchiKey(mol),
            "canonical_smiles": Chem.MolToSmiles(mol),
            "mw": round(Descriptors.MolWt(mol), 3),
            "heavy_atoms": mol.GetNumHeavyAtoms(),
            "rotatable_bonds": rdMolDescriptors.CalcNumRotatableBonds(mol),
        },
        "provenance": record("rdkit", __import__("rdkit").__version__, {"op": "inspect"}),
    }


def project_create(params, srv):
    p = Project.create(_req(params, "path"), params.get("name", "Untitled"))
    return {"result": p.manifest(), "artifacts": [str(p.manifest_path)]}


def project_import_input(params, srv):
    p = Project.open(_path(params))
    return {"result": p.import_input(_path(params, "source"))}


def _lig_opts(params):
    from .ligprep import LigandPrepOptions
    try:
        return LigandPrepOptions(**(params.get("options") or {}))
    except TypeError as e:
        raise InvalidInput(f"unknown ligand option: {e}") from e


def ligand_prepare(params, srv):
    from rdkit import Chem
    from .ligprep import prepare_ligand
    from . import formats
    opts = _lig_opts(params)
    if "path" in params:
        mols, warns = formats.read_molecules(_path(params))
        results = [r for i, m in enumerate(mols)
                   for r in prepare_ligand(m, m.GetProp("_Name") if m.HasProp("_Name") and m.GetProp("_Name") else f"ligand_{i + 1}", opts)]
    elif "smiles" in params:
        warns, results = [], prepare_ligand(_text(params, "smiles"), str(params.get("name", "ligand")), opts)
    else:
        raise InvalidInput("smiles or path required")
    artifacts = []
    outdir = params.get("outdir")
    if outdir:
        from pathlib import Path
        Path(outdir).mkdir(parents=True, exist_ok=True)
        for r in results:
            if r.status != "failed":
                (Path(outdir) / f"{r.name}.pdbqt").write_text(r.pdbqt)
                w = Chem.SDWriter(str(Path(outdir) / f"{r.name}.sdf")); w.write(r.mol); w.close()
                artifacts += [str(Path(outdir) / f"{r.name}.pdbqt"), str(Path(outdir) / f"{r.name}.sdf")]
    return {"result": [r.summary() for r in results], "warnings": warns, "artifacts": artifacts,
            "provenance": record("dockchem.ligprep", __version__, params.get("options") or {})}


def library_prepare(params, srv):
    from pathlib import Path
    from .ligprep import prepare_library, read_records
    from rdkit import Chem
    lib = prepare_library(read_records(_path(params)), _lig_opts(params))
    artifacts = []
    if params.get("outdir"):
        out = Path(_req(params, "outdir")); out.mkdir(parents=True, exist_ok=True)
        for p in lib["prepared"]:
            (out / f"{p.name}.pdbqt").write_text(p.pdbqt); artifacts.append(str(out / f"{p.name}.pdbqt"))
        import json
        (out / "library_report.json").write_text(json.dumps(lib["rows"], indent=2))
        artifacts.append(str(out / "library_report.json"))
    return {"result": {"rows": lib["rows"], "counts": lib["counts"]}, "artifacts": artifacts,
            "provenance": record("dockchem.ligprep", __version__, params.get("options") or {})}


def receptor_prepare(params, srv):
    import json
    from pathlib import Path
    from .recprep import ReceptorPrepOptions, prepare_receptor
    try:
        opts = ReceptorPrepOptions(**(params.get("options") or {}))
    except TypeError as e:
        raise InvalidInput(f"unknown receptor option: {e}") from e
    path, outdir = _path(params), _req(params, "outdir")           # validate before announcing any progress
    if srv:
        srv.progress(params.get("job_id"), 0.0, "preparing receptor")
    res = prepare_receptor(path, outdir, opts)
    report = Path(_req(params, "outdir")) / "prep_report.json"
    report.write_text(json.dumps(res.report, indent=2, default=str))
    arts = [str(res.pdb), str(res.pdbqt), str(report)] + ([str(res.flex_pdbqt)] if res.flex_pdbqt else [])
    arts += [l["path"] for l in res.ligands]
    return {"result": {"ligands": res.ligands, "report": res.report}, "warnings": res.report["warnings"],
            "artifacts": arts, "provenance": record("dockchem.recprep", __version__, params.get("options") or {})}


def ligand_from_structure(params, srv):
    from pathlib import Path
    from rdkit import Chem
    from . import ligand_templates as lt
    r = lt.ligand_from_structure(_path(params), _req(params, "resname"), params.get("chain"), params.get("seq"),
                                 params.get("smiles"), params.get("charge"))
    artifacts = []
    if params.get("outdir"):
        out = Path(_req(params, "outdir")); out.mkdir(parents=True, exist_ok=True)
        sdf = out / f"{params['resname']}.sdf"
        w = Chem.SDWriter(str(sdf)); w.write(r.mol); w.close()
        artifacts.append(str(sdf))
    return {"result": {"source": r.source, "smiles": Chem.MolToSmiles(Chem.RemoveHs(r.mol)),
                       "heavy_atoms": r.mol.GetNumHeavyAtoms(), "audit": r.audit},
            "warnings": r.warnings, "artifacts": artifacts,
            "provenance": record("dockchem.ligand_templates", __version__, {"resname": _req(params, "resname"),
                                                                            "source": r.source})}


def _session(params):
    import gemmi
    from .selection import EditSession
    try:
        st = gemmi.read_structure(_path(params))
    except Exception as e:  # noqa: BLE001
        raise InvalidInput(f"cannot read structure: {e}") from e
    s = EditSession(st, int(params.get("model", 1)))
    s.apply_all(list(params.get("edits") or []))
    return s


def structure_assess(params, srv):
    from .assess import assess
    a = assess(_path(params), int(params.get("model", 1)))
    return {"result": a, "artifacts": [], "provenance": record("dockchem.assess", __version__, {"model": params.get("model", 1)})}


def selection_validate(params, srv):
    from .selection import _split, parse
    cmd = params.get("command")
    if cmd is not None:
        _split(cmd)
        parse(_split(cmd)[1])
    else:
        parse(_req(params, "expr"))
    return {"result": {"ok": True}}


def selection_preview(params, srv):
    """Stateless: the app sends the edit list so far (its undo stack) plus the command being considered."""
    s = _session(params)
    e = s.preview(_req(params, "command"))
    return {"result": {"atoms": e.n_atoms, "residues": e.n_residues, "residue_ids": e.residues,
                       "truncated": e.n_residues > len(e.residues), "warnings": e.warnings,
                       "current": s.summary()}}


def selection_apply(params, srv):
    from . import formats
    s = _session(params)
    out = s.materialize()
    artifacts = []
    if params.get("out"):
        formats.write_structure(out, params["out"])
        artifacts.append(str(params["out"]))
    return {"result": {"summary": s.summary(), "edits": [
        {"command": e.command, "atoms": e.n_atoms, "residues": e.n_residues, "warnings": e.warnings}
        for e in s.edits]}, "artifacts": artifacts,
            "provenance": record("dockchem.selection", __version__, {"edits": s.commands()})}


def quickprep_run(params, srv):
    from .quickprep import quick_prep
    path, outdir = _path(params, "receptor"), _req(params, "outdir")
    ligand = _req(params, "ligand")
    if isinstance(ligand, str) and len(ligand) < 4096 and "\n" not in ligand:
        import os
        if os.path.sep in ligand and not os.path.exists(ligand):
            raise InvalidInput(f"ligand file not found: {ligand}")
    res = quick_prep(path, ligand, outdir, params.get("preset", "standard"), params.get("site"), int(params.get("seed", 42)),
                     bool(params.get("allow_unusable", False)), params.get("receptor_options"), params.get("ligand_options"))
    return {"result": res, "warnings": res["warnings"], "artifacts": [res["engine"]["config"], res["receptor"]["pdbqt"]] + res["ligands"]["files"],
            "provenance": record("dockchem.quickprep", __version__, {"preset": res["preset"], "site": res["site"]["how"]})}


def complex_image2d(params, srv):
    """Profile + 2D interaction image for one docked pose (receptor: prepared PDB with hydrogens)."""
    import re
    from pathlib import Path
    from rdkit import Chem
    from . import formats, interactions, render2d
    rec, pose, out = _path(params, "receptor"), _path(params, "pose"), Path(_req(params, "out"))
    idx = int(params.get("pose_index", 1))
    style = params.get("style", "reference")
    if style not in render2d.THEMES:
        raise InvalidInput(f"style must be one of {', '.join(render2d.THEMES)}")
    affinity = params.get("affinity")
    if pose.lower().endswith(".pdbqt"):
        mols, _ = formats.read_pdbqt_poses(pose)
        scores = [float(x) for x in re.findall(r"REMARK VINA RESULT:\s+(-?\d+\.\d+)", Path(pose).read_text())]
        if affinity is None and len(scores) >= idx:
            affinity = scores[idx - 1]
    else:
        mols = [m for m in Chem.SDMolSupplier(pose, removeHs=False) if m is not None]
    if not 1 <= idx <= len(mols):
        raise InvalidInput(f"pose_index {idx} is out of range (the file has {len(mols)} poses)")
    prof = interactions.profile(Path(rec).read_text(), mols[idx - 1], params.get("target_name", ""),
                                params.get("compound", "ligand"), affinity, params.get("score_label", "AutoDock Vina binding affinity"))
    out.parent.mkdir(parents=True, exist_ok=True)
    n = render2d.draw(prof, out, style=style, dpi=int(params.get("dpi", 300)))
    js = out.with_suffix(".json")
    js.write_text(json.dumps(prof, indent=1))
    return {"result": {"counts": interactions.counts(prof), "residues_drawn": n, "style": style},
            "artifacts": [str(out), str(js)],
            "provenance": record("dockchem.render2d", __version__, {"style": style, "pose_index": idx})}


def format_convert(params, srv):
    from . import formats
    return formats.convert(_path(params, "source"), _req(params, "dest"), params.get("charge", 0))


def format_sniff(params, srv):
    from . import formats
    return {"result": formats.sniff(_path(params))}


def poses_to_sdf(params, srv):
    from . import formats
    r = formats.poses_to_sdf(_path(params, "pdbqt"), _req(params, "sdf"))
    return {"result": {"poses": r["poses"]}, "artifacts": r["artifacts"]}


METHODS = {
    "system.ping": ping,
    "system.versions": versions,
    "test.sleep": sleep_job,
    "ligand.inspect": ligand_inspect,
    "ligand.prepare": ligand_prepare,
    "library.prepare": library_prepare,
    "receptor.prepare": receptor_prepare,
    "ligand.from_structure": ligand_from_structure,
    "system.resources": lambda params, srv: {"result": __import__("dockchem.resources", fromlist=["x"]).describe()},
    "structure.assess": structure_assess,
    "quickprep.run": quickprep_run,
    "complex.image2d": complex_image2d,
    "selection.validate": selection_validate,
    "selection.preview": selection_preview,
    "selection.apply": selection_apply,
    "format.convert": format_convert,
    "format.sniff": format_sniff,
    "poses.to_sdf": poses_to_sdf,
    "project.create": project_create,
    "project.import_input": project_import_input,
}
