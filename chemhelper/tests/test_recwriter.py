"""The Meeko-free receptor PDBQT writer must agree with Meeko atom for atom, and refuse to guess."""
import gemmi
import pytest

from dockchem import recprep, recwriter
from dockchem.errors import InvalidInput, ToolFailure
from dockchem.recprep import ReceptorPrepOptions, prepare_receptor

from fixtures import decorated_structure, peptide_structure


def typed(path):
    out = {}
    for l in open(path):
        if l.startswith("ATOM"):
            name = l[12:16].strip()
            out[(l[17:20].strip(), l[21], l[22:27].strip(), "H" if name == "H1" else name)] = l[77:79].strip()
    return out


@pytest.fixture(scope="module")
def pair(tmp_path_factory):
    src = tmp_path_factory.mktemp("in") / "d.pdb"
    decorated_structure().write_pdb(str(src))
    m = prepare_receptor(src, tmp_path_factory.mktemp("m"), ReceptorPrepOptions(reference_ligand="LIG", pdbqt_writer="meeko"))
    t = prepare_receptor(src, tmp_path_factory.mktemp("t"), ReceptorPrepOptions(reference_ligand="LIG", pdbqt_writer="table"))
    return m, t


def test_table_writer_matches_meeko_atom_for_atom(pair):
    m, t = pair
    a, b = typed(m.pdbqt), typed(t.pdbqt)
    assert a and set(a) == set(b)
    assert {k: v for k, v in a.items() if b[k] != v} == {}


def test_table_run_is_recorded_in_report(pair):
    _, t = pair
    assert any(x["tool"] == "dockchem.recwriter" for x in t.report["tools"])
    assert not any(x["tool"] == "meeko" for x in t.report["tools"])
    assert t.report["status"] == "ok"


def test_table_writer_refuses_flexible_residues(tmp_path):
    p = tmp_path / "p.pdb"
    peptide_structure().write_pdb(str(p))
    with pytest.raises(InvalidInput, match="flexible"):
        prepare_receptor(p, tmp_path / "o", ReceptorPrepOptions(pdbqt_writer="table", flex_residues=["A:3"]))
    with pytest.raises(InvalidInput, match="pdbqt_writer"):
        prepare_receptor(p, tmp_path / "o2", ReceptorPrepOptions(pdbqt_writer="x"))


def test_auto_falls_back_to_table_writer_when_meeko_fails(tmp_path, monkeypatch):
    def boom(*a, **k):
        raise ToolFailure("Meeko receptor preparation failed")
    monkeypatch.setattr(recprep, "_meeko", boom)
    p = tmp_path / "p.pdb"
    peptide_structure().write_pdb(str(p))
    res = prepare_receptor(p, tmp_path / "o", ReceptorPrepOptions())
    assert res.pdbqt.exists() and typed(res.pdbqt)
    assert any("built-in table writer" in w for w in res.report["warnings"])
    with pytest.raises(ToolFailure):
        prepare_receptor(p, tmp_path / "o2", ReceptorPrepOptions(pdbqt_writer="meeko"))
    with pytest.raises(ToolFailure):          # flexible residues cannot fall back
        prepare_receptor(p, tmp_path / "o3", ReceptorPrepOptions(flex_residues=["A:3"]))


def pdb_line(serial, name, res, num, xyz, el, chain="A"):
    return (f"ATOM  {serial:5d} {name:<4s} {res:>3s} {chain}{num:4d}    {xyz[0]:8.3f}{xyz[1]:8.3f}{xyz[2]:8.3f}"
            f"  1.00  0.00          {el:>2s}\n")


def test_unknown_atom_drops_residue_by_name_never_guesses(tmp_path):
    text = (pdb_line(1, "N", "ALA", 1, (0, 0, 0), "N") + pdb_line(2, "CA", "ALA", 1, (1.4, 0, 0), "C")
            + pdb_line(3, "QX", "ALA", 1, (2, 1, 0), "C") + pdb_line(4, "N", "GLY", 2, (5, 0, 0), "N"))
    r = recwriter.write_pdbqt(text, tmp_path / "o.pdbqt")
    assert r["dropped"] == [{"residue": "A:ALA1", "reason": "no atom-type rule for ALA.QX"}]
    assert [l[17:20] for l in open(tmp_path / "o.pdbqt")] == ["GLY"]


def test_unknown_residue_dropped(tmp_path):
    r = recwriter.write_pdbqt(pdb_line(1, "C1", "XYZ", 7, (0, 0, 0), "C"), tmp_path / "o.pdbqt")
    assert r["dropped"][0]["reason"] == "no atom-type rules for XYZ"


def cyx(serial0, num, sg):
    return (pdb_line(serial0, "CB", "CYX", num, (sg[0], sg[1] + 1.8, sg[2]), "C")
            + pdb_line(serial0 + 1, "SG", "CYX", num, sg, "S"))


@pytest.mark.parametrize("sep,want", [(2.03, "SA"), (3.5, "S")])
def test_cyx_sulfur_type_depends_on_disulfide(tmp_path, sep, want):
    text = cyx(1, 1, (0, 0, 0)) + cyx(3, 2, (sep, 0, 0))
    recwriter.write_pdbqt(text, tmp_path / "o.pdbqt")
    assert {l[77:79].strip() for l in open(tmp_path / "o.pdbqt") if l[12:16].strip() == "SG"} == {want}


def test_hydrogen_rule_and_terminal_oxygen(tmp_path):
    text = (pdb_line(1, "N", "GLY", 1, (0, 0, 0), "N") + pdb_line(2, "H1", "GLY", 1, (0.5, 0.9, 0), "H")   # on N
            + pdb_line(3, "CA", "GLY", 1, (1.4, 0, 0), "C") + pdb_line(4, "HA2", "GLY", 1, (1.6, -0.6, 0.8), "H")  # on C
            + pdb_line(5, "OXT", "GLY", 1, (2.5, 1.0, 0), "O"))
    r = recwriter.write_pdbqt(text, tmp_path / "o.pdbqt")
    got = {l[12:16].strip(): l[77:79].strip() for l in open(tmp_path / "o.pdbqt")}
    assert got == {"N": "N", "H1": "HD", "CA": "C", "OXT": "OA"} and r["merged_hydrogens"] == 1


def test_output_format_is_valid_pdbqt_columns(pair):
    _, t = pair
    for l in open(t.pdbqt):
        assert l.startswith("ATOM") and len(l.rstrip("\n")) == 79
        float(l[30:38]); float(l[38:46]); float(l[46:54]); float(l[70:76])
