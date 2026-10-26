"""The `cores` setting: default 2, validated, clamped, and applied everywhere the helper runs in parallel."""
import os

import pytest

from dockchem import resources
from dockchem.errors import InvalidInput
from dockchem.server import Server


def rpc(method, **params):
    return Server().handle({"jsonrpc": "2.0", "id": 1, "method": method, "params": params})


def test_default_is_two_cores(monkeypatch):
    monkeypatch.delenv("DOCKCHEM_CORES", raising=False)
    assert resources.cores() == min(2, resources.total_cores())
    assert resources.DEFAULT_CORES == 2


def test_environment_variable_sets_the_default(monkeypatch):
    monkeypatch.setenv("DOCKCHEM_CORES", "1")
    assert resources.cores() == 1
    monkeypatch.setenv("DOCKCHEM_CORES", "garbage")
    assert resources.cores() == min(2, resources.total_cores())


@pytest.mark.parametrize("bad", [0, -1, 1.5, "two", None, True, [2], {"n": 2}, "", "0", float("nan"), float("inf"), "-inf"])
def test_invalid_core_counts_are_rejected(bad):
    with pytest.raises(InvalidInput, match="whole number"):
        resources.validate(bad)


def test_more_cores_than_the_machine_has_are_clamped_with_a_warning():
    n, warn = resources.validate(resources.total_cores() + 50)
    assert n == resources.total_cores() and "using" in warn
    assert resources.validate(1) == (1, None) and resources.validate("1") == (1, None)


def test_a_request_can_set_cores_and_gets_it_back_in_provenance():
    r = rpc("system.ping", cores=1)
    assert r["result"]["provenance"]["cores"] == 1 and "error" not in r
    assert rpc("system.ping")["result"]["provenance"]["cores"] == resources.cores()


def test_bad_cores_in_a_request_is_invalid_input_and_clamping_warns():
    assert rpc("system.ping", cores=0)["error"]["data"]["type"] == "InvalidInput"
    assert rpc("system.ping", cores="many")["error"]["data"]["type"] == "InvalidInput"
    r = rpc("system.ping", cores=resources.total_cores() + 9)
    assert any("using" in w for w in r["result"]["warnings"])


def test_the_setting_does_not_leak_between_requests():
    rpc("system.ping", cores=1)
    assert resources.cores() == min(2, resources.total_cores()) or os.environ.get("DOCKCHEM_CORES")


def test_resources_report_for_the_ui_slider():
    d = rpc("system.resources")["result"]["result"]
    assert d["min"] == 1 and d["max"] == resources.total_cores() and d["default"] == min(2, d["max"])
    assert 1 <= d["performance_cores"] <= d["total_cores"] and "current" in d


def test_helper_subprocesses_get_the_thread_limit(tmp_path):
    import sys
    from dockchem import recprep
    h = resources.bind(1)
    try:
        r = recprep._run([sys.executable, "-c", "import os;print(os.environ['OMP_NUM_THREADS'], os.environ['OPENBLAS_NUM_THREADS'])"], tmp_path)
    finally:
        resources.unbind(h)
    assert r.stdout.split() == ["1", "1"]


def test_rdkit_embedding_uses_the_requested_threads(monkeypatch):
    from rdkit.Chem import AllChem
    from dockchem.ligprep import LigandPrepOptions, prepare_ligand
    seen = []
    real = AllChem.EmbedMultipleConfs

    def spy(mol, n, ps):
        seen.append(ps.numThreads)
        return real(mol, n, ps)
    monkeypatch.setattr(AllChem, "EmbedMultipleConfs", spy)
    h = resources.bind(1)
    try:
        prepare_ligand("CCOC(=O)c1ccccc1", "x", LigandPrepOptions())
    finally:
        resources.unbind(h)
    assert seen and set(seen) == {1}


def test_openmm_minimisation_runs_with_the_requested_threads(tmp_path, monkeypatch):
    import openmm as mm
    from openmm import app
    from dockchem import recprep
    from fixtures import peptide_structure
    from dockchem.recprep import ReceptorPrepOptions, prepare_receptor
    src = tmp_path / "p.pdb"
    peptide_structure().write_pdb(str(src))
    base = prepare_receptor(src, tmp_path / "b", ReceptorPrepOptions(crosscheck=False))
    props = []
    real = app.Simulation

    def spy(top, system, integ, platform=None, platformProperties=None):
        props.append(platformProperties)
        return real(top, system, integ, platform, platformProperties)
    monkeypatch.setattr(app, "Simulation", spy)
    h = resources.bind(1)
    try:
        recprep._minimize(base.pdb, tmp_path / "o.pdb", recprep._Report())
    finally:
        resources.unbind(h)
    assert props == [{"Threads": "1"}]


def test_quick_prep_writes_the_core_count_into_the_vina_config(tmp_path):
    from dockchem.quickprep import quick_prep
    from fixtures import decorated_structure
    src = tmp_path / "d.pdb"
    decorated_structure().write_pdb(str(src))
    r = quick_prep(src, "CCO", tmp_path / "default")
    assert f"cpu = {resources.cores()}" in open(r["engine"]["config"]).read() and r["engine"]["cores"] == resources.cores()
    h = resources.bind(1)
    try:
        r1 = quick_prep(src, "CCO", tmp_path / "one")
    finally:
        resources.unbind(h)
    assert "cpu = 1" in open(r1["engine"]["config"]).read() and r1["engine"]["cores"] == 1
