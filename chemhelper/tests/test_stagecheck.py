"""Per-stage output checks: each must fail loudly on a bad result and stay quiet on legitimate changes."""
import pytest

from dockchem import stagecheck as S
from dockchem.errors import ToolFailure


def pdb(serial, name, res, num, xyz, el, chain="A"):
    return (f"ATOM  {serial:5d} {name:<4s} {res:>3s} {chain}{num:4d}    {xyz[0]:8.3f}{xyz[1]:8.3f}{xyz[2]:8.3f}"
            f"  1.00  0.00          {el:>2s}\n")


def pdbqt(serial, name, res, num, xyz, ad, chain="A"):
    return (f"ATOM  {serial:5d} {name:<4s} {res:>3s} {chain}{num:4d}    {xyz[0]:8.3f}{xyz[1]:8.3f}{xyz[2]:8.3f}"
            f"  1.00  0.00    +0.000 {ad:<2s}\n")


BEFORE = pdb(1, "N", "ALA", 1, (0, 0, 0), "N") + pdb(2, "CA", "ALA", 1, (1.4, 0, 0), "C") + pdb(3, "N", "GLY", 2, (5, 0, 0), "N")


def test_unchanged_atoms_pass_and_report_zero_movement():
    assert S.check_moved_atoms("pdb2pqr", BEFORE, BEFORE) == 0.0


def test_lost_residue_is_a_stage_failure():
    after = pdb(1, "N", "ALA", 1, (0, 0, 0), "N") + pdb(2, "CA", "ALA", 1, (1.4, 0, 0), "C")
    with pytest.raises(ToolFailure, match="residues lost: A:2") as e:
        S.check_moved_atoms("pdb2pqr", BEFORE, after)
    assert e.value.stage == "pdb2pqr"


def test_atom_moved_too_far_fails_unless_exempt():
    after = BEFORE.replace("   1.400", "   3.000")
    with pytest.raises(ToolFailure, match="A:ALA1.CA moved 1.60"):
        S.check_moved_atoms("pdb2pqr", BEFORE, after)
    assert S.check_moved_atoms("pdb2pqr", BEFORE, after, allowed_moved={("A", "1")}) == 0.0     # exempt: not counted
    S.check_moved_atoms("pdb2pqr", BEFORE, BEFORE.replace("   1.400", "   1.900"))     # 0.5 A: within the limit


def test_hydrogens_moving_is_not_an_issue():
    b = BEFORE + pdb(4, "HA", "ALA", 1, (1.5, 1.0, 0), "H")
    a = BEFORE + pdb(4, "HA", "ALA", 1, (9.5, 1.0, 0), "H")
    S.check_moved_atoms("pdb2pqr", b, a)


def test_protonated_model_needs_hydrogens_on_every_amino_acid():
    with pytest.raises(ToolFailure, match="A:ALA1 has no hydrogens"):
        S.check_hydrogens("pdb2pqr", BEFORE)
    S.check_hydrogens("pdb2pqr", BEFORE + pdb(4, "H", "ALA", 1, (0, 1, 0), "H") + pdb(5, "H", "GLY", 2, (5, 1, 0), "H"))


def prepared():
    return BEFORE.replace("N   GLY", "N   GLY")


def test_pdbqt_complete_passes(tmp_path):
    q = tmp_path / "r.pdbqt"
    q.write_text(pdbqt(1, "N", "ALA", 1, (0, 0, 0), "N") + pdbqt(2, "CA", "ALA", 1, (1.4, 0, 0), "C")
                 + pdbqt(3, "N", "GLY", 2, (5, 0, 0), "N"))
    assert S.check_pdbqt(BEFORE, q, []) == 3


def test_pdbqt_missing_residue_must_be_reported(tmp_path):
    q = tmp_path / "r.pdbqt"
    q.write_text(pdbqt(1, "N", "ALA", 1, (0, 0, 0), "N") + pdbqt(2, "CA", "ALA", 1, (1.4, 0, 0), "C"))
    with pytest.raises(ToolFailure, match="GLY2 is absent from the PDBQT without being reported"):
        S.check_pdbqt(BEFORE, q, [])
    assert S.check_pdbqt(BEFORE, q, [{"residue": "A:GLY2", "reason": "x"}]) == 2
    assert S.check_pdbqt(BEFORE, q, [{"residue": "A:2", "reason": "Meeko-style key"}]) == 2


def test_pdbqt_lost_atoms_unknown_types_and_empty(tmp_path):
    q = tmp_path / "r.pdbqt"
    q.write_text(pdbqt(1, "N", "ALA", 1, (0, 0, 0), "N") + pdbqt(3, "N", "GLY", 2, (5, 0, 0), "N"))
    with pytest.raises(ToolFailure, match="lost heavy atoms"):
        S.check_pdbqt(BEFORE, q, [])
    q.write_text(pdbqt(1, "N", "ALA", 1, (0, 0, 0), "QQ") + pdbqt(2, "CA", "ALA", 1, (1.4, 0, 0), "C")
                 + pdbqt(3, "N", "GLY", 2, (5, 0, 0), "N"))
    with pytest.raises(ToolFailure, match="unknown atom type 'QQ'"):
        S.check_pdbqt(BEFORE, q, [])
    q.write_text("")
    with pytest.raises(ToolFailure):
        S.check_pdbqt("", q, [])


def test_pdbqt_tolerates_renamed_water_atoms_oxt_and_flex_file(tmp_path):
    prep = BEFORE + pdb(4, "O", "HOH", 9, (8, 8, 8), "O") + pdb(5, "OXT", "GLY", 2, (6, 1, 0), "O")
    q, f = tmp_path / "r.pdbqt", tmp_path / "f.pdbqt"
    q.write_text(pdbqt(1, "N", "ALA", 1, (0, 0, 0), "N") + pdbqt(3, "N", "GLY", 2, (5, 0, 0), "N")
                 + pdbqt(4, "OW", "HOH", 9, (8, 8, 8), "OA"))            # water atom renamed; OXT omitted
    f.write_text(pdbqt(2, "CA", "ALA", 1, (1.4, 0, 0), "C"))              # flexible part lives elsewhere
    assert S.check_pdbqt(prep, q, [], f) == 4
