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


def test_backbone_anchor_moved_fails_unless_exempt():
    after = BEFORE.replace("   1.400", "   3.000")
    with pytest.raises(ToolFailure, match="A:ALA1.CA moved 1.60"):
        S.check_moved_atoms("pdb2pqr", BEFORE, after)
    assert S.check_moved_atoms("pdb2pqr", BEFORE, after, allowed_moved={("A", "1")}) == 0.0     # exempt: not counted
    S.check_moved_atoms("pdb2pqr", BEFORE, BEFORE.replace("   1.400", "   1.700"), limit=0.5)  # 0.3 A: fine


def residue(res, atoms, num=1):
    return "".join(pdb(i + 1, n, res, num, xyz, n[0]) for i, (n, xyz) in enumerate(atoms))


GLU = [("N", (0, 0, 0)), ("CA", (1.45, 0, 0)), ("C", (2.0, 1.4, 0)), ("CB", (2.0, -0.8, 1.2)), ("CG", (3.5, -0.8, 1.2)),
       ("CD", (4.1, -2.1, 1.2)), ("OE1", (5.3, -2.2, 1.2)), ("OE2", (3.4, -3.2, 1.2))]


def test_equivalent_atom_relabelling_and_side_chain_rotation_are_legitimate():
    """PDB2PQR swaps OE1/OE2 and rotates side chains to relieve clashes; neither changes any bond length."""
    swapped = [(n if n not in ("OE1", "OE2") else {"OE1": "OE2", "OE2": "OE1"}[n], xyz) for n, xyz in GLU]
    S.check_moved_atoms("pdb2pqr", residue("GLU", GLU), residue("GLU", swapped))
    rotated = [(n, (x, -y - 1.6 if n in ("CG", "CD", "OE1", "OE2") else y, z)) for n, (x, y, z) in GLU]   # mirror about CB
    rotated = [(n, (x, 2 * -0.8 - y, z) if n in ("CG", "CD", "OE1", "OE2") else (x, y, z)) for n, (x, y, z) in GLU]
    S.check_moved_atoms("pdb2pqr", residue("GLU", GLU), residue("GLU", rotated))


def test_a_stretched_or_squashed_bond_is_a_stage_failure():
    broken = [(n, (x + 1.2, y, z) if n == "CD" else (x, y, z)) for n, (x, y, z) in GLU]
    with pytest.raises(ToolFailure, match="bond CD-OE1 changed|bond CG-CD changed"):
        S.check_moved_atoms("pdb2pqr", residue("GLU", GLU), residue("GLU", broken))


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


def test_missing_amide_hydrogen_is_listed_but_not_for_proline_or_when_present():
    t = (pdb(1, "N", "LYS", 74, (0, 0, 0), "N") + pdb(2, "CA", "LYS", 74, (1.4, 0, 0), "C")             # no H at all
         + pdb(3, "N", "ALA", 75, (5, 0, 0), "N") + pdb(4, "H", "ALA", 75, (5, 1, 0), "H")                # fine
         + pdb(5, "N", "PRO", 76, (9, 0, 0), "N")                                                          # Pro has none
         + pdb(6, "N", "GLY", 77, (13, 0, 0), "N") + pdb(7, "H1", "GLY", 77, (13, 1, 0), "H"))           # ammonium H
    assert S.missing_amide_hydrogens(t) == ["A:LYS74"]
