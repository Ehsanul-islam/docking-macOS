"""Per-residue protonation overrides and flip accept/reject (gap 4)."""
import math
from pathlib import Path

import gemmi
import pytest

from dockchem import protonation as P
from dockchem.errors import InvalidInput, ToolFailure
from dockchem.recprep import ReceptorPrepOptions, prepare_receptor

from fixtures import peptide_structure

REPO = Path(__file__).resolve().parents[2]
PDB_3PTB = REPO / "pdb_structures_50" / "3PTB.pdb"


def atom(serial, name, res, num, xyz, chain="A"):
    return (f"ATOM  {serial:5d} {name:<4s} {res:>3s} {chain}{num:4d}    "
            f"{xyz[0]:8.3f}{xyz[1]:8.3f}{xyz[2]:8.3f}  1.00 20.00\n")


ASN = [("CB", (0, 0, 0)), ("CG", (1.5, 0, 0)), ("OD1", (2.2, 1.1, 0)), ("ND2", (2.2, -1.1, 0)),
       ("HD21", (3.1, -1.2, 0.0)), ("HD22", (1.7, -1.9, 0.0))]


def asn_text(swap=False, res="ASN"):
    rows = []
    for i, (n, xyz) in enumerate(ASN):
        x, y, z = xyz
        if swap and n not in ("CB", "CG"):
            y = -y                      # mirror through the CB-CG axis == half turn about it (z = 0 here)
        rows.append(atom(i + 1, n, res, 25, (x, y, z)))
    return "".join(rows)


def pos(text, name):
    for l in text.splitlines():
        if l[12:16].strip() == name:
            return (float(l[30:38]), float(l[38:46]), float(l[46:54]))


# ---- key parsing and validation -------------------------------------------------------------------

@pytest.mark.parametrize("bad", ["45", "A", "AB:4", ":4", "A:", "A:x", "A:-"])
def test_bad_keys_rejected(bad):
    with pytest.raises(InvalidInput):
        P.parse_key(bad)


def test_key_with_insertion_code_and_negative_number():
    assert P.parse_key("A:45B") == ("A", "45B")
    assert P.parse_key("B:-3") == ("B", "-3")


def test_override_unknown_residue():
    with pytest.raises(InvalidInput, match="no such residue"):
        P.apply_overrides(atom(1, "CA", "HIS", 3, (0, 0, 0)), {"A:9": "HIE"})


def test_override_state_must_belong_to_residue():
    with pytest.raises(InvalidInput, match="not a state of HIS"):
        P.apply_overrides(atom(1, "CA", "HIS", 3, (0, 0, 0)), {"A:3": "ASH"})


def test_override_on_residue_without_states():
    with pytest.raises(InvalidInput, match="no selectable"):
        P.apply_overrides(atom(1, "CA", "GLY", 3, (0, 0, 0)), {"A:3": "HIE"})


def test_cyx_cannot_be_chosen():
    with pytest.raises(InvalidInput):
        P.apply_overrides(atom(1, "CA", "CYS", 3, (0, 0, 0)), {"A:3": "CYX"})


def test_override_is_case_insensitive_and_renames_every_atom():
    text = atom(1, "N", "HIS", 3, (0, 0, 0)) + atom(2, "CA", "HIS", 3, (1, 0, 0))
    out, done = P.apply_overrides(text, {"A:3": "hie"})
    assert done == [{"key": "A:3", "old": "HIS", "new": "HIE"}]
    assert all(l[17:20] == "HIE" for l in out.splitlines())


def test_verify_detects_override_that_did_not_take_effect():
    applied = [{"key": "A:3", "old": "HIS", "new": "HIE"}]
    wrong_name = atom(1, "N", "HID", 3, (0, 0, 0))
    with pytest.raises(ToolFailure, match="not honoured"):
        P.verify_overrides(wrong_name, applied)
    wrong_h = atom(1, "N", "HIE", 3, (0, 0, 0)) + atom(2, "HD1", "HIE", 3, (1, 0, 0))
    with pytest.raises(ToolFailure, match="hydrogens do not match"):
        P.verify_overrides(wrong_h, applied)
    right = atom(1, "N", "HIE", 3, (0, 0, 0)) + atom(2, "HE2", "HIE", 3, (1, 0, 0))
    P.verify_overrides(right, applied)


# ---- flips ----------------------------------------------------------------------------------------

def test_flip_detected_when_terminal_atoms_swap():
    assert P.detect_flips(asn_text(), asn_text(swap=True)) == [{"key": "A:25", "resname": "ASN"}]


def test_no_flip_when_unchanged_or_only_moved_slightly():
    assert P.detect_flips(asn_text(), asn_text()) == []
    assert P.detect_flips(asn_text(), asn_text().replace("2.200", "2.250")) == []


def test_undo_flip_restores_original_coordinates_including_hydrogens():
    original = asn_text()
    restored = P.undo_flips(asn_text(swap=True), ["A:25"])
    for name, xyz in ASN:
        assert math.dist(pos(restored, name), xyz) < 0.01, name


def test_undo_flip_rejects_non_flippable_and_unknown_residues():
    with pytest.raises(InvalidInput, match="not an Asn"):
        P.undo_flips(atom(1, "CA", "ALA", 5, (0, 0, 0)), ["A:5"])
    with pytest.raises(InvalidInput, match="no such residue"):
        P.undo_flips(asn_text(), ["A:99"])


# ---- end to end -----------------------------------------------------------------------------------

@pytest.fixture(scope="module")
def peptide(tmp_path_factory):
    p = tmp_path_factory.mktemp("pep") / "p.pdb"
    peptide_structure().write_pdb(str(p))
    return p


def final_names(res):
    st = gemmi.read_pdb(str(res.pdb))
    return {r.seqid.num: r.name for r in st[0][0]}


def test_overrides_reach_the_prepared_receptor_and_report(peptide, tmp_path):
    want = {"A:3": "HIE", "A:4": "ASH", "A:5": "GLH", "A:6": "LYN"}
    res = prepare_receptor(peptide, tmp_path, ReceptorPrepOptions(protonation_overrides=want))
    names = final_names(res)
    assert [names[3], names[4], names[5], names[6]] == ["HIE", "ASH", "GLH", "LYN"]
    ch = {c["residue"]: c for c in res.report["changes"] if c["stage"] == "protonation"}
    assert ch["A:HIS3"]["reason"] == "user override" and ch["A:HIS3"]["tool"] == "user"
    assert any("set by user override" in w for w in res.report["warnings"])


def test_default_run_is_unchanged_without_overrides(peptide, tmp_path):
    res = prepare_receptor(peptide, tmp_path, ReceptorPrepOptions())
    assert final_names(res)[3] == "HID"
    assert not any(c["reason"] == "user override" for c in res.report["changes"])


def test_override_on_missing_residue_is_an_error_not_ignored(peptide, tmp_path):
    with pytest.raises(InvalidInput, match="no such residue"):
        prepare_receptor(peptide, tmp_path, ReceptorPrepOptions(protonation_overrides={"A:99": "HIE"}))


def test_override_on_residue_removed_by_edit_is_an_error(peptide, tmp_path):
    with pytest.raises(InvalidInput, match="no such residue"):
        prepare_receptor(peptide, tmp_path, ReceptorPrepOptions(
            edits=["delete resi 3"], protonation_overrides={"A:3": "HIE"}))


def test_bad_flip_decision_value(peptide, tmp_path):
    with pytest.raises(InvalidInput, match="accept"):
        prepare_receptor(peptide, tmp_path, ReceptorPrepOptions(flip_decisions={"A:9": "maybe"}))


def test_rejecting_a_flip_that_did_not_happen_warns(peptide, tmp_path):
    res = prepare_receptor(peptide, tmp_path, ReceptorPrepOptions(flip_decisions={"A:9": "reject"}))
    assert any("nothing to reject" in w for w in res.report["warnings"])


@pytest.mark.skipif(not PDB_3PTB.exists(), reason="real structure not available")
def test_real_structure_flips_are_listed_and_can_be_rejected(tmp_path):
    def heavy(path, num, names):
        st = gemmi.read_pdb(str(path))
        for r in st[0]["A"]:
            if r.seqid.num == num and r.name == "ASN":
                return {n: gemmi.Position(r[n][0].pos.x, r[n][0].pos.y, r[n][0].pos.z) for n in names}   # copy: the structure dies here

    src = heavy(PDB_3PTB, 25, ("OD1", "ND2"))
    accepted = prepare_receptor(PDB_3PTB, tmp_path / "a", ReceptorPrepOptions(reference_ligand="BEN"))
    flips = [c for c in accepted.report["changes"] if c["stage"] == "flip"]
    assert "A:ASN25" in {c["residue"] for c in flips}
    a = heavy(accepted.pdb, 25, ("OD1", "ND2"))
    assert a["OD1"].dist(src["ND2"]) < 0.5            # accepted: O and N traded places

    rejected = prepare_receptor(PDB_3PTB, tmp_path / "r", ReceptorPrepOptions(
        reference_ligand="BEN", flip_decisions={"A:25": "reject"}))
    r = heavy(rejected.pdb, 25, ("OD1", "ND2"))
    assert r["OD1"].dist(src["OD1"]) < 0.5 and r["ND2"].dist(src["ND2"]) < 0.5
    rc = {c["residue"]: c for c in rejected.report["changes"] if c["stage"] == "flip"}
    assert rc["A:ASN25"]["tool"] == "user" and "restored" in rc["A:ASN25"]["change"]
    assert rc["A:GLN30"]["change"] == "side chain flipped"          # other flips untouched
