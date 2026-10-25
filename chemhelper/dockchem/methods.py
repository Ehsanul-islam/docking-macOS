"""RPC method registry. Phase 0 exposes system + project + a ligand smoke method."""
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
    if not smi:
        raise InvalidInput("smiles required")
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
    p = Project.create(params["path"], params.get("name", "Untitled"))
    return {"result": p.manifest(), "artifacts": [str(p.manifest_path)]}


def project_import_input(params, srv):
    p = Project.open(params["path"])
    return {"result": p.import_input(params["source"])}


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
        mols, warns = formats.read_molecules(params["path"])
        results = [r for i, m in enumerate(mols)
                   for r in prepare_ligand(m, m.GetProp("_Name") if m.HasProp("_Name") and m.GetProp("_Name") else f"ligand_{i + 1}", opts)]
    elif "smiles" in params:
        warns, results = [], prepare_ligand(params["smiles"], params.get("name", "ligand"), opts)
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
    lib = prepare_library(read_records(params["path"]), _lig_opts(params))
    artifacts = []
    if params.get("outdir"):
        out = Path(params["outdir"]); out.mkdir(parents=True, exist_ok=True)
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
    if srv:
        srv.progress(params.get("job_id"), 0.0, "preparing receptor")
    res = prepare_receptor(params["path"], params["outdir"], opts)
    report = Path(params["outdir"]) / "prep_report.json"
    report.write_text(json.dumps(res.report, indent=2, default=str))
    arts = [str(res.pdb), str(res.pdbqt), str(report)] + ([str(res.flex_pdbqt)] if res.flex_pdbqt else [])
    arts += [l["path"] for l in res.ligands]
    return {"result": {"ligands": res.ligands, "report": res.report}, "warnings": res.report["warnings"],
            "artifacts": arts, "provenance": record("dockchem.recprep", __version__, params.get("options") or {})}


def ligand_from_structure(params, srv):
    from pathlib import Path
    from rdkit import Chem
    from . import ligand_templates as lt
    r = lt.ligand_from_structure(params["path"], params["resname"], params.get("chain"), params.get("seq"),
                                 params.get("smiles"), params.get("charge"))
    artifacts = []
    if params.get("outdir"):
        out = Path(params["outdir"]); out.mkdir(parents=True, exist_ok=True)
        sdf = out / f"{params['resname']}.sdf"
        w = Chem.SDWriter(str(sdf)); w.write(r.mol); w.close()
        artifacts.append(str(sdf))
    return {"result": {"source": r.source, "smiles": Chem.MolToSmiles(Chem.RemoveHs(r.mol)),
                       "heavy_atoms": r.mol.GetNumHeavyAtoms(), "audit": r.audit},
            "warnings": r.warnings, "artifacts": artifacts,
            "provenance": record("dockchem.ligand_templates", __version__, {"resname": params["resname"],
                                                                            "source": r.source})}


def _session(params):
    import gemmi
    from .selection import EditSession
    try:
        st = gemmi.read_structure(str(params["path"]))
    except Exception as e:  # noqa: BLE001
        raise InvalidInput(f"cannot read structure: {e}") from e
    s = EditSession(st, int(params.get("model", 1)))
    s.apply_all(list(params.get("edits") or []))
    return s


def structure_assess(params, srv):
    from .assess import assess
    a = assess(params["path"], int(params.get("model", 1)))
    return {"result": a, "artifacts": [], "provenance": record("dockchem.assess", __version__, {"model": params.get("model", 1)})}


def selection_validate(params, srv):
    from .selection import _split, parse
    cmd = params.get("command")
    if cmd is not None:
        _split(cmd)
        parse(_split(cmd)[1])
    else:
        parse(params["expr"])
    return {"result": {"ok": True}}


def selection_preview(params, srv):
    """Stateless: the app sends the edit list so far (its undo stack) plus the command being considered."""
    s = _session(params)
    e = s.preview(params["command"])
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


def format_convert(params, srv):
    from . import formats
    return formats.convert(params["source"], params["dest"], params.get("charge", 0))


def format_sniff(params, srv):
    from . import formats
    return {"result": formats.sniff(params["path"])}


def poses_to_sdf(params, srv):
    from . import formats
    r = formats.poses_to_sdf(params["pdbqt"], params["sdf"])
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
    "structure.assess": structure_assess,
    "selection.validate": selection_validate,
    "selection.preview": selection_preview,
    "selection.apply": selection_apply,
    "format.convert": format_convert,
    "format.sniff": format_sniff,
    "poses.to_sdf": poses_to_sdf,
    "project.create": project_create,
    "project.import_input": project_import_input,
}
