"""Batch docking runner: scheduling, ranking, failures, timeout, cancel, resume, validation."""
import io
import json
import os
import stat
import threading
import time
from pathlib import Path

import pytest

from dockchem import cancel, docking, resources
from dockchem.errors import Cancelled, InvalidInput, UnsupportedFeature

FAKE = '''#!/usr/bin/env python3
import sys, time, pathlib
a = sys.argv[1:]
if "--version" in a:
    print("AutoDock Vina vFAKE"); sys.exit(0)
g = lambda k: a[a.index(k) + 1]
lig, out, cpu = g("--ligand"), g("--out"), g("--cpu")
name = pathlib.Path(lig).stem
with open(pathlib.Path(out).parent.parent / "invocations.txt", "a") as f:
    f.write(f"{name} {cpu} {g('--exhaustiveness')}\\n")
if "fail" in name:
    print("Parse error on line 3"); sys.exit(1)
if "slow" in name:
    time.sleep(60)
if "nopose" in name:
    sys.exit(0)
aff = -{"a": 9.1, "b": 7.5, "c": 8.3, "d": 6.0}.get(name[0], 5.5)
open(out, "w").write(f"MODEL 1\\nREMARK VINA RESULT:    {aff:.3f}      0.000      0.000\\nATOM      1  C   LIG A   1       0.000   0.000   0.000  1.00  0.00     0.000 C\\nENDMDL\\n"
                     f"MODEL 2\\nREMARK VINA RESULT:    {aff + 1:.3f}      1.200      2.500\\nATOM      1  C   LIG A   1       0.000   0.000   0.000  1.00  0.00     0.000 C\\nENDMDL\\n")
print("done")
'''
LIG = "ROOT\nATOM      1  C   LIG A   1       0.000   0.000   0.000  1.00  0.00     0.000 C\nENDROOT\nTORSDOF 0\n"


@pytest.fixture
def engine(tmp_path, monkeypatch):
    e = tmp_path / "fake_vina"
    e.write_text(FAKE)
    e.chmod(e.stat().st_mode | stat.S_IEXEC)
    monkeypatch.setenv("DOCKCHEM_VINA", str(e))
    return e


@pytest.fixture
def setup(tmp_path, engine):
    rec = tmp_path / "rec.pdbqt"
    rec.write_text("ATOM      1  C   ALA A   1       0.000   0.000   0.000  1.00  0.00     0.000 C\n")
    def make(*names):
        out = []
        for n in names:
            f = tmp_path / "ligs" / f"{n}.pdbqt"
            f.parent.mkdir(exist_ok=True)
            f.write_text(LIG)
            out.append(str(f))
        return out
    return rec, make


def run(rec, ligs, out, **kw):
    kw.setdefault("write_sdf", False)
    return docking.run_batch(rec, ligs, out, [0, 0, 0], [20, 20, 20], **kw)


def invocations(out):
    p = Path(out) / "invocations.txt"
    return [l.split() for l in p.read_text().splitlines()] if p.exists() else []


def test_results_are_ranked_best_first_and_written_everywhere(tmp_path, setup):
    rec, make = setup
    r = run(rec, make("a_lig", "b_lig", "c_lig", "d_lig"), tmp_path / "o", exhaustiveness=4, num_modes=2)
    assert [x["name"] for x in r["results"]] == ["a_lig", "c_lig", "b_lig", "d_lig"]
    assert [x["rank"] for x in r["results"]] == [1, 2, 3, 4] and r["results"][0]["affinity"] == -9.1
    assert r["docked"] == 4 and r["failed"] == 0 and "not a measured free energy" in r["scores_are"]
    rows = (tmp_path / "o" / "results.csv").read_text().splitlines()
    assert rows[0].startswith("rank,name,affinity_kcal_mol") and rows[1].startswith("1,a_lig,-9.1")
    run_json = json.loads((tmp_path / "o" / "run.json").read_text())
    assert run_json["settings"]["engine"] == "AutoDock Vina vFAKE" and run_json["settings"]["exhaustiveness"] == 4
    assert (tmp_path / "o" / "poses" / "a_lig_out.pdbqt").exists() and (tmp_path / "o" / "logs" / "a_lig.log").exists()


@pytest.mark.parametrize("n,cores,expect", [(100, 4, (4, 1)), (4, 4, (4, 1)), (1, 4, (1, 4)), (2, 4, (2, 2)), (3, 8, (3, 2)), (5, 1, (1, 1))])
def test_core_budget_is_split_between_processes(n, cores, expect):
    assert docking.plan_cores(n, cores) == expect


def test_the_engine_gets_the_planned_threads(tmp_path, setup):
    rec, make = setup
    h = resources.bind(min(4, resources.total_cores()))
    try:
        run(rec, make("a1", "b1"), tmp_path / "two")
        run(rec, make("a2"), tmp_path / "one")
    finally:
        resources.unbind(h)
    c = min(4, resources.total_cores())
    assert {i[1] for i in invocations(tmp_path / "one")} == {str(c)}
    assert {i[1] for i in invocations(tmp_path / "two")} == {str(docking.plan_cores(2, c)[1])}


def test_a_failing_or_silent_ligand_never_stops_the_batch(tmp_path, setup):
    rec, make = setup
    r = run(rec, make("a_ok", "fail_x", "nopose_y", "b_ok"), tmp_path / "o")
    by = {x["name"]: x for x in r["results"]}
    assert by["a_ok"]["status"] == by["b_ok"]["status"] == "docked"
    assert by["fail_x"]["status"] == "failed" and "Parse error" in by["fail_x"]["error"] and "code 1" in by["fail_x"]["error"]
    assert by["nopose_y"]["status"] == "failed" and "no poses" in by["nopose_y"]["error"]
    assert r["docked"] == 2 and r["failed"] == 2
    assert [x["name"] for x in r["results"]][:2] == ["a_ok", "b_ok"]                   # failures listed after the ranking


def test_a_ligand_over_the_time_limit_is_stopped_and_reported(tmp_path, setup):
    rec, make = setup
    t0 = time.time()
    r = run(rec, make("slow_one", "a_ok"), tmp_path / "o", timeout_s=1)
    assert time.time() - t0 < 15
    by = {x["name"]: x for x in r["results"]}
    assert by["slow_one"]["status"] == "failed" and "exceeded the 1 s limit" in by["slow_one"]["error"] and by["a_ok"]["status"] == "docked"


def test_cancel_kills_running_engines_and_keeps_finished_results(tmp_path, setup):
    rec, make = setup
    tok = cancel.Token()
    h = cancel.bind(tok)
    resources_h = resources.bind(1)                    # one at a time: a_first finishes, slow_second is running when cancelled
    try:
        threading.Timer(1.5, tok.cancel).start()
        t0 = time.time()
        with pytest.raises(Cancelled, match="after 1 of 3 ligands"):
            run(rec, make("a_first", "slow_second", "b_third"), tmp_path / "o")
        assert time.time() - t0 < 10
    finally:
        resources.unbind(resources_h)
        cancel.unbind(h)
    assert "a_first" in (tmp_path / "o" / "results.csv").read_text()
    time.sleep(0.3)
    me = str(os.getpid())
    alive = [p for p in os.listdir("/proc") if p.isdigit() and p != me and
             b"slow_second" in open(f"/proc/{p}/cmdline", "rb").read() if os.access(f"/proc/{p}/cmdline", os.R_OK)]
    assert alive == [], [open(f"/proc/{p}/cmdline", "rb").read()[:200] for p in alive]


def test_a_rerun_resumes_and_a_changed_setting_is_refused(tmp_path, setup):
    rec, make = setup
    ligs = make("a_lig", "b_lig")
    run(rec, ligs, tmp_path / "o", exhaustiveness=8)
    before = len(invocations(tmp_path / "o"))
    again = run(rec, make("a_lig", "b_lig", "c_lig"), tmp_path / "o", exhaustiveness=8)
    assert len(invocations(tmp_path / "o")) == before + 1 and again["docked"] == 3          # only the new ligand ran
    with pytest.raises(InvalidInput, match="different settings"):
        run(rec, ligs, tmp_path / "o", exhaustiveness=16)
    fresh = run(rec, ligs, tmp_path / "o", exhaustiveness=16, resume=False)
    assert fresh["docked"] == 2


def test_same_file_name_in_different_folders_does_not_collide(tmp_path, setup):
    rec, make = setup
    one = make("a_same")[0]
    other = tmp_path / "other" / "a_same.pdbqt"
    other.parent.mkdir()
    other.write_text(LIG)
    r = run(rec, [one, str(other)], tmp_path / "o")
    assert sorted(x["name"] for x in r["results"]) == ["a_same", "a_same_2"]


def test_bad_inputs_are_refused_before_any_engine_runs(tmp_path, setup, monkeypatch):
    rec, make = setup
    ok = make("a_ok")
    for kwargs, msg in [(dict(ligands=[]), "no ligands"), (dict(ligands=[str(tmp_path / "nope.pdbqt")]), "not found"),
                        (dict(ligands=ok, size=[0, 5, 5]), "positive"), (dict(ligands=ok, exhaustiveness=0), "exhaustiveness"),
                        (dict(ligands=ok, num_modes=99), "num_modes"), (dict(ligands=ok, scoring="ad4"), "not supported")]:
        args = dict(receptor=rec, ligands=ok, outdir=tmp_path / "x", center=[0, 0, 0], size=[20, 20, 20])
        args.update(kwargs)
        with pytest.raises(InvalidInput, match=msg):
            docking.run_batch(**args)
    bad = tmp_path / "bad.pdbqt"
    bad.write_text("ATOM garbage\n")
    with pytest.raises(InvalidInput, match="not a ligand PDBQT"):
        run(rec, [str(bad)], tmp_path / "y")
    with pytest.raises(InvalidInput, match="receptor PDBQT not found"):
        run(tmp_path / "missing.pdbqt", ok, tmp_path / "z")
    assert not invocations(tmp_path / "x")


def test_a_missing_engine_is_a_clear_unsupported_feature(monkeypatch):
    monkeypatch.delenv("DOCKCHEM_VINA", raising=False)
    monkeypatch.setenv("PATH", "/nonexistent")
    monkeypatch.setattr(docking, "__file__", "/nonexistent/a/b/c.py")
    with pytest.raises(UnsupportedFeature, match="DOCKCHEM_VINA"):
        docking.find_engine()


def test_rpc_runs_a_quick_prep_folder_and_reports_progress(tmp_path, engine):
    from dockchem.quickprep import quick_prep
    from dockchem.server import Server
    from fixtures import decorated_structure
    src = tmp_path / "d.pdb"
    decorated_structure().write_pdb(str(src))
    lib = tmp_path / "lib.smi"
    lib.write_text("CCO a_ethanol\nCCCO b_propanol\n")
    quick_prep(src, str(lib), tmp_path / "qp", preset="fast")
    buf = io.StringIO()
    resp = Server(stdout=buf).handle({"jsonrpc": "2.0", "id": 1, "method": "dock.run",
                                      "params": {"quickprep_dir": str(tmp_path / "qp"), "outdir": str(tmp_path / "run"),
                                                 "job_id": "J", "write_sdf": False}})
    res = resp["result"]["result"]
    assert res["docked"] == 2 and res["results"][0]["rank"] == 1
    prog = [json.loads(l) for l in buf.getvalue().splitlines()]
    assert [p["params"]["fraction"] for p in prog] == [0.5, 1.0] and "kcal/mol" in prog[0]["params"]["message"]
    assert resp["result"]["provenance"]["tool"] == "dockchem.docking"
    bad = Server().handle({"jsonrpc": "2.0", "id": 2, "method": "dock.run", "params": {"quickprep_dir": str(tmp_path), "outdir": str(tmp_path / "r2")}})
    assert bad["error"]["data"]["type"] == "InvalidInput"


REAL = os.environ.get("DOCKCHEM_VINA") or "/tmp/claude-0/vina/vina_1.2.7_linux_x86_64"


@pytest.mark.skipif(not Path(REAL).is_file(), reason="real Vina binary not available")
def test_real_vina_end_to_end_on_a_tiny_search(tmp_path, monkeypatch):
    from dockchem.quickprep import quick_prep
    from fixtures import decorated_structure
    monkeypatch.setenv("DOCKCHEM_VINA", REAL)
    src = tmp_path / "d.pdb"
    decorated_structure().write_pdb(str(src))
    lib = tmp_path / "lib.smi"
    lib.write_text("CCO ethanol\nCC(=O)N acetamide\n")
    q = quick_prep(src, str(lib), tmp_path / "qp", preset="fast")
    r = docking.run_batch(q["receptor"]["pdbqt"], q["ligands"]["files"], tmp_path / "run", q["site"]["center"], q["site"]["size"],
                          exhaustiveness=1, num_modes=2)
    assert r["docked"] == 2 and all(x["affinity"] < 0 and x["n_poses"] >= 1 for x in r["results"])
    assert all(Path(x["sdf"]).exists() for x in r["results"])
