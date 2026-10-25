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


METHODS = {
    "system.ping": ping,
    "system.versions": versions,
    "test.sleep": sleep_job,
    "ligand.inspect": ligand_inspect,
    "project.create": project_create,
    "project.import_input": project_import_input,
}
