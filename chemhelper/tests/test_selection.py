"""Selection language and non-destructive editing. Real structures from pdb_structures_50/ plus oracles that
do not use the code under test (direct gemmi loops, brute-force distances)."""
import json
import subprocess
import sys
import time
from pathlib import Path

import gemmi
import numpy as np
import pytest

from dockchem import selection as S
from dockchem.errors import InvalidInput, SelectionError

PDBS = Path(__file__).resolve().parents[2] / "pdb_structures_50"
pytestmark = pytest.mark.skipif(not PDBS.is_dir(), reason="pdb_structures_50 not present")


@pytest.fixture(scope="module")
def iep():
    st = gemmi.read_structure(str(PDBS / "1IEP.pdb"))
    return st, S.build_table(st)


def count(t, expr):
    return int(S.select(t, expr).sum())


def atoms(st):
    return [(ch.name, r.name, r.seqid.num, a.name, a.element.name, a.b_iso, a.pos)
            for ch in st[0] for r in ch for a in r]


# ---------------------------------------------------------------- terms against direct oracles
def test_table_matches_structure(iep):
    st, t = iep
    assert t.n == sum(len(r) for ch in st[0] for r in ch) == 4710


@pytest.mark.parametrize("expr,oracle", [
    ("all", lambda a: True),
    ("chain A", lambda a: a[0] == "A"),
    ("chain A,B", lambda a: a[0] in "AB"),
    ("resn HOH", lambda a: a[1] == "HOH"),
    ("resn STI", lambda a: a[1] == "STI"),
    ("name CA", lambda a: a[3] == "CA"),
    ("elem S", lambda a: a[4] == "S"),
    ("elem S,N", lambda a: a[4] in ("S", "N")),
    ("resi 300-310", lambda a: 300 <= a[2] <= 310),
    ("resi 300-310 and name CA", lambda a: 300 <= a[2] <= 310 and a[3] == "CA"),
    ("chain A and not resn HOH", lambda a: a[0] == "A" and a[1] != "HOH"),
    ("bfactor > 60", lambda a: a[5] > 60),
    ("name C*", lambda a: a[3].startswith("C")),
    ("resn S?I", lambda a: len(a[1]) == 3 and a[1][0] == "S" and a[1][2] == "I"),
])
def test_atom_level_terms(iep, expr, oracle):
    st, t = iep
    assert count(t, expr) == sum(1 for a in atoms(st) if oracle(a)), expr


def test_class_words(iep):
    _, t = iep
    assert count(t, "water") == 172 and count(t, "resn STI") == count(t, "ligand") == 74
    assert count(t, "hetero") == count(t, "water") + count(t, "ligand") + count(t, "additive") + count(t, "metal")
    assert count(t, "protein") + count(t, "hetero") == t.n
    assert count(t, "backbone") + count(t, "sidechain") + count(t, "protein and hydrogen") \
        == count(t, "protein")


def test_within_and_around_match_brute_force(iep):
    st, t = iep
    lig = t.xyz[S.select(t, "resn STI")]
    d = np.linalg.norm(t.xyz[:, None, :] - lig[None, :, :], axis=2).min(axis=1)
    assert count(t, "within 5 of resn STI") == int((d <= 5).sum())
    assert count(t, "around 5 of resn STI") == int(((d <= 5) & ~S.select(t, "resn STI")).sum())
    assert count(t, "within 3 of resn STI") <= count(t, "within 5 of resn STI")


def test_byres_selects_whole_residues(iep):
    _, t = iep
    sel = S.select(t, "byres (within 4 of resn STI) and protein")
    for g in np.unique(t.res_id[sel]):
        s, e = t.res_slices[g]
        assert sel[s:e].all()
    assert sel.sum() > count(t, "within 4 of resn STI and protein")


# ---------------------------------------------------------------- grammar
def test_precedence_and_parentheses(iep):
    _, t = iep
    a, b, w = (S.select(t, x) for x in ("chain A", "chain B", "water"))
    assert (S.select(t, "chain A or chain B and water") == (a | (b & w))).all()
    assert (S.select(t, "(chain A or chain B) and water") == ((a | b) & w)).all()
    assert (S.select(t, "not chain A and water") == (~a & w)).all()
    assert (S.select(t, "not (chain A and water)") == ~(a & w)).all()


def test_keywords_are_case_insensitive_and_whitespace_tolerant(iep):
    _, t = iep
    assert count(t, "  CHAIN   A  AND   WATER ") == count(t, "chain A and water")
    assert count(t, "(water)") == count(t, "water")


@pytest.mark.parametrize("expr,needle", [
    ("", "empty"), ("chain", "needs a value"), ("chain and water", "needs a value"),
    ("(water", "missing ')'"), ("water)", "unexpected"), ("water and", "ends unexpectedly"),
    ("flurb A", "unknown keyword"), ("within 5 resn STI", "expected 'of'"), ("within x of water", "distance"),
    ("within 500 of water", "between 0 and 50"), ("bfactor 5", "comparison"), ("bfactor > x", "number"),
    ("resi 5-1", "backwards"), ("resi abc", "bad residue range"),
])
def test_bad_expressions_give_located_errors(iep, expr, needle):
    _, t = iep
    with pytest.raises(SelectionError) as e:
        S.select(t, expr)
    assert needle in e.value.message and e.value.expr == expr and 0 <= e.value.pos <= len(expr)


# ---------------------------------------------------------------- edit session
def test_delete_keep_and_counts(iep):
    st, _ = iep
    s = S.EditSession(st)
    e = s.apply("delete water")
    assert (e.n_atoms, e.n_residues) == (172, 172) and e.warnings == []
    s.apply("keep chain A")
    out = s.materialize()
    assert [ch.name for ch in out[0]] == ["A"]
    assert all(r.name != "HOH" for r in out[0]["A"])
    assert s.summary()["atoms"] == sum(len(r) for ch in out[0] for r in ch)


def test_input_structure_is_never_modified(iep):
    st, _ = iep
    before = [(a[0], a[2], a[3]) for a in atoms(st)]
    s = S.EditSession(st)
    s.apply_all(["delete chain B", "delete water"])
    s.materialize()
    assert [(a[0], a[2], a[3]) for a in atoms(st)] == before


def test_preview_does_not_change_state(iep):
    st, _ = iep
    s = S.EditSession(st)
    s.preview("delete protein")
    assert s.summary()["atoms"] == 4710 and s.edits == []


def test_undo_redo_by_replay_is_exact(iep):
    st, _ = iep
    cmds = ["delete water", "delete chain B", "delete resn STI"]
    full = S.EditSession(st); full.apply_all(cmds)
    undone = S.EditSession(st); undone.apply_all(cmds[:-1])                  # undo = drop the last command
    redone = S.EditSession(st); redone.apply_all(cmds[:-1] + cmds[-1:])       # redo = put it back
    assert undone.summary()["atoms"] == full.summary()["atoms"] + 37
    assert (redone.alive == full.alive).all()


def test_replay_is_deterministic(iep, tmp_path):
    st, _ = iep
    outs = []
    for i in range(2):
        s = S.EditSession(st); s.apply_all(["delete chain B", "delete within 3 of resn STI and water"])
        p = tmp_path / f"o{i}.pdb"; s.materialize().write_pdb(str(p)); outs.append(p.read_text())
    assert outs[0] == outs[1]


def test_materialized_structure_roundtrips_through_files(iep, tmp_path):
    st, _ = iep
    s = S.EditSession(st); s.apply_all(["delete water", "delete chain B"])
    out = s.materialize()
    for ext in ("pdb", "cif"):
        p = tmp_path / f"o.{ext}"
        (out.write_pdb if ext == "pdb" else out.make_mmcif_document().write_file)(str(p))
        back = gemmi.read_structure(str(p))
        assert sum(len(r) for ch in back[0] for r in ch) == s.summary()["atoms"]


def test_partial_residue_deletion_is_warned_about(iep):
    st, _ = iep
    e = S.EditSession(st).preview("delete name CA")                          # atoms, not whole residues
    assert any("lose only some atoms" in w for w in e.warnings)


def test_whole_residue_deletion_has_no_partial_warning(iep):
    st, _ = iep
    e = S.EditSession(st).preview("delete water")
    assert not any("lose only some" in w for w in e.warnings)


def test_pocket_and_chain_break_and_majority_warnings(iep):
    st, _ = iep
    e = S.EditSession(st).preview("delete byres (within 4 of resn STI) and protein and chain A")
    assert any("within 5 A of a ligand" in w for w in e.warnings)
    e2 = S.EditSession(st).preview("delete byres name CA and resi 300-305 and chain A")
    assert any("new chain break" in w for w in e2.warnings)
    e3 = S.EditSession(st).preview("delete protein")
    assert any("every protein residue" in w for w in e3.warnings)
    e4 = S.EditSession(st).preview("delete protein and not resi 300-310")      # well over half
    assert any("protein residues" in w and " of " in w for w in e4.warnings)
    e5 = S.EditSession(st).preview("delete chain A and protein")                # exactly half: no majority warning
    assert not any("protein residues" in w for w in e5.warnings)


def test_bad_edit_verb(iep):
    st, _ = iep
    with pytest.raises(SelectionError):
        S.EditSession(st).apply("zap water")


def test_multi_model_nmr_structure_uses_requested_model():
    st = gemmi.read_structure(str(PDBS / "1L2Y.pdb"))
    s = S.EditSession(st, model=3)
    s.apply("delete hydrogen")
    out = s.materialize()
    assert len(out) == 1
    heavy = [a for ch in st[2] for r in ch for a in r if a.element.name != "H"]
    got = [a for ch in out[0] for r in ch for a in r]
    assert len(got) == len(heavy) and abs(got[0].pos.x - heavy[0].pos.x) < 1e-6
    with pytest.raises(InvalidInput):
        S.EditSession(st, model=99)


def test_large_structure_is_fast():
    st = gemmi.read_structure(str(PDBS / "1IGT.pdb"))            # ~13k atoms, 4 chains
    t0 = time.time()
    s = S.EditSession(st)
    s.apply_all(["delete water", "delete within 4 of resn FUL and water", "keep chain A,B"])
    s.materialize()
    assert time.time() - t0 < 10


# ---------------------------------------------------------------- RPC and receptor prep
def rpc(*reqs):
    p = subprocess.run([sys.executable, "-m", "dockchem.server"], capture_output=True, text=True, timeout=120,
                       input="\n".join(json.dumps({"jsonrpc": "2.0", "id": i, "method": m, "params": pr})
                                       for i, (m, pr) in enumerate(reqs, 1)) + "\n")
    return {m["id"]: m for m in map(json.loads, p.stdout.splitlines())}


def test_rpc_validate_preview_apply_and_errors(tmp_path):
    path = str(PDBS / "1IEP.pdb")
    out = rpc(("selection.validate", {"expr": "chain A and water"}),
              ("selection.validate", {"expr": "chain A and"}),
              ("selection.preview", {"path": path, "edits": ["delete water"], "command": "delete chain B"}),
              ("selection.apply", {"path": path, "edits": ["delete water", "delete chain B"],
                                   "out": str(tmp_path / "o.pdb")}))
    assert out[1]["result"]["result"]["ok"] is True
    err = out[2]["error"]["data"]
    assert err["type"] == "SelectionError" and err["pos"] >= 0 and err["expr"] == "chain A and"
    pv = out[3]["result"]["result"]
    assert pv["current"]["atoms"] == 4710 - 172 and pv["atoms"] == 2268 and pv["residues"] == 277
    ap = out[4]["result"]
    assert ap["result"]["summary"]["chains"] == ["A"] and (tmp_path / "o.pdb").exists()
    assert ap["provenance"]["params"]["edits"] == ["delete water", "delete chain B"]


def test_receptor_prep_applies_edits_and_logs_them(tmp_path):
    from dockchem.recprep import ReceptorPrepOptions, prepare_receptor
    r = prepare_receptor(PDBS / "1STP.pdb", tmp_path / "o", ReceptorPrepOptions(
        reference_ligand="BTN", edits=["delete resn SO4,GOL", "delete within 8 of resn BTN and water"]))
    log = [c for c in r.report["changes"] if c["stage"] == "selection"]
    assert [c["change"] for c in log] == ["delete resn SO4,GOL", "delete within 8 of resn BTN and water"]
    assert r.ligands and r.ligands[0]["resname"] == "BTN" and r.pdbqt.exists()


def test_receptor_prep_edits_can_drop_a_chain(tmp_path):
    from dockchem.recprep import ReceptorPrepOptions, prepare_receptor
    r = prepare_receptor(PDBS / "1AKI.pdb", tmp_path / "o", ReceptorPrepOptions(
        edits=["delete water", "delete resi 1-10"]))
    seqs = {int(l[22:26]) for l in r.pdb.read_text().splitlines() if l.startswith("ATOM")}
    assert min(seqs) == 11 and 5 not in seqs


def test_receptor_prep_refuses_edits_that_remove_everything_or_the_reference(tmp_path):
    from dockchem.recprep import ReceptorPrepOptions, prepare_receptor
    with pytest.raises(InvalidInput, match="every protein residue"):
        prepare_receptor(PDBS / "1AKI.pdb", tmp_path / "a", ReceptorPrepOptions(edits=["delete protein"]))
    with pytest.raises(InvalidInput, match="removed by one of the edits"):
        prepare_receptor(PDBS / "1STP.pdb", tmp_path / "b",
                         ReceptorPrepOptions(reference_ligand="BTN", edits=["delete resn BTN"]))
    with pytest.raises(SelectionError):
        prepare_receptor(PDBS / "1AKI.pdb", tmp_path / "c", ReceptorPrepOptions(edits=["delete chain"]))
