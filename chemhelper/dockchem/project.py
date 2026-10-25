""".dockproj package store (plan §4.6, §14)."""
import json
import shutil
import uuid
from pathlib import Path

from .errors import InvalidInput
from .provenance import record, sha256_file

SCHEMA_VERSION = 1
SUBDIRS = ("inputs", "receptor", "ligands", "runs", "reports")


class Project:
    def __init__(self, root: Path):
        self.root = Path(root)

    @property
    def manifest_path(self) -> Path:
        return self.root / "project.json"

    @property
    def provenance_path(self) -> Path:
        return self.root / "provenance.jsonl"

    @classmethod
    def create(cls, root: str | Path, name: str) -> "Project":
        root = Path(root)
        if root.suffix != ".dockproj":
            raise InvalidInput("project path must end in .dockproj")
        if root.exists():
            raise InvalidInput(f"{root} already exists")
        for d in SUBDIRS:
            (root / d).mkdir(parents=True)
        p = cls(root)
        p._write_manifest({
            "schema_version": SCHEMA_VERSION,
            "id": str(uuid.uuid4()),
            "name": name,
            "inputs": [],
            "parameters": {},
            "seeds": {},
        })
        p.provenance_path.touch()
        return p

    @classmethod
    def open(cls, root: str | Path) -> "Project":
        p = cls(root)
        if not p.manifest_path.is_file():
            raise InvalidInput(f"{root} is not a .dockproj package")
        m = p.manifest()
        if m.get("schema_version", 0) > SCHEMA_VERSION:
            raise InvalidInput("project was created by a newer app version")
        return p

    def manifest(self) -> dict:
        return json.loads(self.manifest_path.read_text())

    def _write_manifest(self, m: dict) -> None:
        tmp = self.manifest_path.with_suffix(".tmp")
        tmp.write_text(json.dumps(m, indent=2, sort_keys=True) + "\n")
        tmp.replace(self.manifest_path)

    def log(self, operation: str, prov: dict, artifacts: list[str] | None = None) -> None:
        """Append-only provenance log."""
        entry = {"operation": operation, **prov, "artifacts": artifacts or []}
        with open(self.provenance_path, "a") as f:
            f.write(json.dumps(entry, sort_keys=True) + "\n")

    def import_input(self, src: str | Path) -> dict:
        """Copy a user file into inputs/ untouched, recording its SHA-256."""
        src = Path(src)
        if not src.is_file():
            raise InvalidInput(f"{src} is not a file")
        dest = self.root / "inputs" / src.name
        if dest.exists():
            raise InvalidInput(f"input {src.name} already imported")
        shutil.copy2(src, dest)
        entry = {"path": f"inputs/{src.name}", "sha256": sha256_file(dest)}
        m = self.manifest()
        m["inputs"].append(entry)
        self._write_manifest(m)
        self.log("import_input", record("dockchem", "import", {"source": src.name}),
                 [entry["path"]])
        return entry

    def verify_inputs(self) -> list[str]:
        """Return paths whose checksum no longer matches the manifest."""
        return [e["path"] for e in self.manifest()["inputs"]
                if sha256_file(self.root / e["path"]) != e["sha256"]]
