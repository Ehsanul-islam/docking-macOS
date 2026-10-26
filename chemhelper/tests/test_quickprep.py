"""One-panel quick prep and the image RPC."""
import copy
import json

import pytest

from dockchem import methods
from dockchem.errors import InvalidInput, UnusableStructure
from dockchem.quickprep import PRESETS, quick_prep
from dockchem.server import Server

from fixtures import decorated_structure, peptide_structure


@pytest.fixture(scope="module")
def src(tmp_path_factory):
    p = tmp_path_factory.mktemp("qp") / "d.pdb"
    decorated_structure().write_pdb(str(p))
    return p


@pytest.fixture(scope="module")
def done(src, tmp_path_factory):
    return quick_prep(src, "CC(=O)Nc1ccc(O)cc1", tmp_path_factory.mktemp("o"), preset="standard")


def test_everything_needed_to_dock_is_produced_and_every_decision_is_reported(done):
    assert done["ready_to_dock"] and done["assessment"]["verdict"] in ("ready", "fixable")
    assert "LIG" in done["site"]["how"] and all(s >= 20 for s in done["site"]["size"])
    assert done["receptor"]["status"] in ("ok", "degraded") and done["ligands"]["prepared"] == 1
    conf = open(done["engine"]["config"]).read()
    assert "exhaustiveness = 32" in conf and "num_modes = 9" in conf and "seed = 42" in conf and "receptor = " in conf
    assert [s["step"] for s in done["steps"]] == ["assess", "site", "receptor", "ligands"]


def test_presets_change_only_the_engine_settings(src, tmp_path):
    r = quick_prep(src, "CCO", tmp_path, preset="thorough")
    assert r["engine"]["exhaustiveness"] == 64 and r["engine"]["num_modes"] == 20
    assert set(PRESETS) == {"fast", "standard", "thorough"}


def test_site_overrides_box_and_residues_and_blind_docking(src, tmp_path):
    c = quick_prep(src, "CCO", tmp_path / "c", site={"center": [1, 2, 3], "size": [18, 19, 20]})
    assert c["site"]["center"] == [1.0, 2.0, 3.0] and c["site"]["size"] == [18.0, 19.0, 20.0] and "given by the user" in c["site"]["how"]
    r = quick_prep(src, "CCO", tmp_path / "r", site={"residues": ["A:5", "A:6"]})
    assert "residues A:5, A:6" in r["site"]["how"]
    peptide = tmp_path / "pep.pdb"
    peptide_structure("GSHDEKRCNQFYTALGSHDEK").write_pdb(str(peptide))
    b = quick_prep(peptide, "CCO", tmp_path / "b")
    assert "blind docking" in b["site"]["how"] and any("blind docking" in w for w in b["warnings"])


def test_bad_inputs_are_refused_clearly(src, tmp_path):
    with pytest.raises(InvalidInput, match="preset must be"):
        quick_prep(src, "CCO", tmp_path, preset="max")
    with pytest.raises(InvalidInput, match="no ligand could be prepared"):
        quick_prep(src, "not a molecule", tmp_path / "x")
    with pytest.raises(InvalidInput, match="SMILES string or the path"):
        quick_prep(src, "", tmp_path / "y")
    with pytest.raises(InvalidInput, match="were found"):
        quick_prep(src, "CCO", tmp_path / "z", site={"residues": ["A:9999"]})


def test_an_unusable_structure_is_refused_with_the_reasons(tmp_path):
    st = copy.deepcopy(peptide_structure("GSHDEKRCNQFYTALGSHDEK"))
    for r in st[0][0]:
        for i in reversed(range(len(r))):
            if r[i].name != "CA":
                del r[i]
    f = tmp_path / "ca.pdb"
    st.write_pdb(str(f))
    with pytest.raises(UnusableStructure, match="Calpha-only"):
        quick_prep(f, "CCO", tmp_path / "o")


def test_a_library_file_is_prepared_molecule_by_molecule(src, tmp_path):
    lib = tmp_path / "lib.smi"
    lib.write_text("CCO ethanol\nc1ccccc1 benzene\nnot_a_smiles bad\n")
    r = quick_prep(src, str(lib), tmp_path / "o")
    assert r["ligands"]["prepared"] == 2 and len(r["ligands"]["failed"]) == 1


def test_rpc_quickprep_and_image_validation(src, tmp_path):
    def call(method, **p):
        return Server().handle({"jsonrpc": "2.0", "id": 1, "method": method, "params": p})
    ok = call("quickprep.run", receptor=str(src), ligand="CCO", outdir=str(tmp_path / "o"), preset="fast")
    assert ok["result"]["result"]["ready_to_dock"] and ok["result"]["artifacts"]
    assert call("quickprep.run", receptor=str(src), outdir=str(tmp_path))["error"]["data"]["type"] == "InvalidInput"
    assert call("complex.image2d", receptor=str(src), pose=str(src), out=str(tmp_path / "x.png"), style="neon")["error"]["data"]["type"] == "InvalidInput"
