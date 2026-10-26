import pytest
import json
import subprocess
import sys


def run(lines):
    p = subprocess.run([sys.executable, "-m", "dockchem.server"],
                       input="\n".join(json.dumps(l) if isinstance(l, dict) else l for l in lines) + "\n",
                       capture_output=True, text=True, timeout=60)
    return [json.loads(l) for l in p.stdout.splitlines()]


def req(i, method, **params):
    return {"jsonrpc": "2.0", "id": i, "method": method, "params": params}


def test_ping_roundtrip():
    out = run([req(1, "system.ping")])
    assert out[0]["id"] == 1
    assert out[0]["result"]["result"] == "pong"
    assert set(out[0]["result"]) >= {"warnings", "provenance", "artifacts"}


def test_versions():
    r = run([req(1, "system.versions")])[0]["result"]["result"]
    assert r["rdkit"] and r["gemmi"]


def test_ligand_inspect():
    r = run([req(1, "ligand.inspect", smiles="CCO")])[0]["result"]
    assert r["result"]["heavy_atoms"] == 3
    assert r["provenance"]["tool"] == "rdkit"


def test_typed_errors():
    out = run([req(1, "ligand.inspect", smiles="not-a-smiles"), req(2, "ligand.inspect")])
    assert out[0]["error"]["data"]["type"] == "InvalidInput"
    assert out[1]["error"]["data"]["type"] == "InvalidInput"


def test_unknown_method_and_parse_error():
    out = run([req(1, "nope"), "{bad json"])
    assert out[0]["error"]["code"] == -32601
    assert out[1]["error"]["code"] == -32700


def test_progress_notifications():
    out = run([req(1, "test.sleep", job_id="j", steps=3)])
    prog = [m for m in out if m.get("method") == "progress"]
    assert [p["params"]["fraction"] for p in prog] == [1 / 3, 2 / 3, 1.0]
    assert out[-1]["result"]["result"] == "done"


def test_server_survives_after_error():
    out = run([req(1, "ligand.inspect", smiles="bad("), req(2, "system.ping")])
    assert "error" in out[0] and out[1]["result"]["result"] == "pong"


def test_prepare_methods_over_rpc(tmp_path):
    import sys; sys.path.insert(0, str(__import__("pathlib").Path(__file__).parent))
    from fixtures import decorated_structure
    rec = tmp_path / "r.pdb"; decorated_structure().write_pdb(str(rec))
    lib = tmp_path / "l.smi"; lib.write_text("CCO eth\nbad bad\n")
    out = run([req(1, "ligand.prepare", smiles="CC(=O)Nc1ccccc1", outdir=str(tmp_path / "lig")),
               req(2, "library.prepare", path=str(lib), outdir=str(tmp_path / "lib")),
               req(3, "receptor.prepare", path=str(rec), outdir=str(tmp_path / "rec"),
                   options={"reference_ligand": "LIG"}),
               req(4, "format.sniff", path=str(rec)),
               req(5, "ligand.prepare", smiles="CCO", options={"nope": 1})])
    by = {m["id"]: m for m in out if "id" in m}
    assert by[1]["result"]["result"][0]["status"] == "prepared"
    assert (tmp_path / "lig" / "ligand.pdbqt").exists()
    assert by[2]["result"]["result"]["counts"]["failed"] == 1
    assert by[3]["result"]["result"]["ligands"][0]["resname"] == "LIG"
    assert (tmp_path / "rec" / "prep_report.json").exists()
    assert by[4]["result"]["result"] == "receptor"
    assert by[5]["error"]["data"]["type"] == "InvalidInput"


# ---- hardening: cancellation, input guards, internal errors ------------------------------------------------

import threading
import time


def test_cancel_during_a_running_job_stops_it_promptly():
    """The reader thread must see `cancel` while the worker is busy (the old single-threaded loop could not)."""
    p = subprocess.Popen([sys.executable, "-m", "dockchem.server"], stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                         stderr=subprocess.PIPE, text=True)
    p.stdin.write(json.dumps(req(1, "test.sleep", job_id="long", steps=200, interval=0.1)) + "\n")
    p.stdin.flush()
    time.sleep(1.0)
    t0 = time.time()
    p.stdin.write(json.dumps({"jsonrpc": "2.0", "method": "cancel", "params": {"job_id": "long"}}) + "\n")
    p.stdin.close()
    out = [json.loads(l) for l in p.stdout.read().splitlines()]
    p.wait(timeout=20)
    final = [m for m in out if m.get("id") == 1][0]
    assert final["error"]["data"]["type"] == "Cancelled"
    assert time.time() - t0 < 5                                       # 200 steps x 0.1 s would have taken 20 s


def test_cancel_before_the_job_starts_cancels_it():
    out = run([{"jsonrpc": "2.0", "method": "cancel", "params": {"job_id": "early"}},
               req(1, "test.sleep", job_id="early", steps=5)])
    assert out[0]["error"]["data"]["type"] == "Cancelled"


def test_responses_stay_in_request_order_with_a_worker_thread():
    out = run([req(1, "system.ping"), "{bad json", req(2, "ligand.inspect", smiles="CCO"), req(3, "nope")])
    assert [m.get("id") for m in out] == [1, None, 2, 3]


def test_missing_and_unreadable_inputs_are_invalid_input_not_internal_errors(tmp_path):
    folder = tmp_path / "d"
    folder.mkdir()
    empty = tmp_path / "empty.pdb"
    empty.write_text("")
    calls = [req(1, "receptor.prepare", outdir=str(tmp_path / "o")),                         # no path at all
             req(2, "receptor.prepare", path=str(tmp_path / "nope.pdb"), outdir=str(tmp_path / "o")),
             req(3, "structure.assess", path=str(folder)),
             req(4, "structure.assess", path=str(empty)),
             req(5, "format.sniff"), req(6, "ligand.inspect", smiles=""), req(7, "selection.validate"),
             req(8, "ligand.prepare", path=str(folder))]
    for r in [m for m in run(calls) if "id" in m]:
        assert r["error"]["data"]["type"] == "InvalidInput", r


def test_a_bug_is_reported_as_internal_error_with_a_traceback(monkeypatch):
    from dockchem import methods
    from dockchem.server import Server
    monkeypatch.setitem(methods.METHODS, "test.bug", lambda params, srv: {}["missing"])
    resp = Server().handle(req(1, "test.bug"))
    assert resp["error"]["data"]["type"] == "InternalError" and "KeyError" in resp["error"]["message"]
    assert "Traceback" in resp["error"]["data"]["log"]


def test_oversized_input_is_rejected_with_the_limit(tmp_path, monkeypatch):
    from dockchem import limits
    monkeypatch.setattr(limits, "MAX_INPUT_BYTES", 10)
    f = tmp_path / "big.pdb"
    f.write_text("ATOM" * 10)
    with pytest.raises(Exception, match="over the .* MB limit"):
        limits.check_input_file(f)
