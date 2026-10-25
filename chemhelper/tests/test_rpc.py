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
