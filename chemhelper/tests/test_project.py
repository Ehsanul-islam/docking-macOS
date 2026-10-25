import json

import pytest

from dockchem.errors import InvalidInput
from dockchem.project import Project


def test_create_layout(tmp_path):
    p = Project.create(tmp_path / "a.dockproj", "demo")
    for d in ("inputs", "receptor", "ligands", "runs", "reports"):
        assert (p.root / d).is_dir()
    assert p.manifest()["schema_version"] == 1
    assert p.provenance_path.exists()


def test_rejects_bad_suffix_and_existing(tmp_path):
    with pytest.raises(InvalidInput):
        Project.create(tmp_path / "a", "x")
    Project.create(tmp_path / "a.dockproj", "x")
    with pytest.raises(InvalidInput):
        Project.create(tmp_path / "a.dockproj", "x")


def test_import_checksum_and_provenance(tmp_path):
    src = tmp_path / "rec.pdb"
    src.write_text("ATOM\n")
    p = Project.create(tmp_path / "a.dockproj", "x")
    e = p.import_input(src)
    assert len(e["sha256"]) == 64
    assert p.verify_inputs() == []
    log = [json.loads(l) for l in p.provenance_path.read_text().splitlines()]
    assert log[0]["operation"] == "import_input" and log[0]["arch"]
    (p.root / e["path"]).write_text("tampered")
    assert p.verify_inputs() == [e["path"]]


def test_newer_schema_refused(tmp_path):
    p = Project.create(tmp_path / "a.dockproj", "x")
    m = p.manifest(); m["schema_version"] = 99
    p.manifest_path.write_text(json.dumps(m))
    with pytest.raises(InvalidInput):
        Project.open(p.root)
