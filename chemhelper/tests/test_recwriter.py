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


def test_auto_uses_the_table_writer_and_is_unaffected_by_meeko_failing(tmp_path, monkeypatch):
    def boom(*a, **k):
        raise ToolFailure("Meeko receptor preparation failed")
    monkeypatch.setattr(recprep, "_meeko", boom)
    p = tmp_path / "p.pdb"
    peptide_structure().write_pdb(str(p))
    res = prepare_receptor(p, tmp_path / "o", ReceptorPrepOptions())
    assert res.pdbqt.exists() and typed(res.pdbqt) and res.report["status"] == "ok"
    assert res.report["pdbqt"]["writer"] == "table"
    assert res.report["pdbqt"]["crosscheck"]["result"] == "unavailable"
    with pytest.raises(ToolFailure):          # explicit Meeko stays strict
        prepare_receptor(p, tmp_path / "o2", ReceptorPrepOptions(pdbqt_writer="meeko"))
    with pytest.raises(ToolFailure):          # flexible residues cannot use the table writer
        prepare_receptor(p, tmp_path / "o3", ReceptorPrepOptions(flex_residues=["A:3"]))


def test_default_run_cross_checks_against_meeko_and_they_agree(tmp_path):
    p = tmp_path / "d.pdb"
    decorated_structure().write_pdb(str(p))
    res = prepare_receptor(p, tmp_path / "o", ReceptorPrepOptions(reference_ligand="LIG"))
    cc = res.report["pdbqt"]["crosscheck"]
    assert res.report["pdbqt"]["writer"] == "table" and cc["result"] == "agree" and cc["n_residues_compared"] > 10
    assert any(t["tool"] == "meeko" and t["params"]["role"] == "cross-check" for t in res.report["tools"])
    assert not (tmp_path / "o" / "meeko_crosscheck").exists()        # scratch is cleaned up


def test_crosscheck_can_be_switched_off(tmp_path):
    p = tmp_path / "p.pdb"
    peptide_structure().write_pdb(str(p))
    res = prepare_receptor(p, tmp_path / "o", ReceptorPrepOptions(crosscheck=False))
    assert "crosscheck" not in res.report["pdbqt"] and not any(t["tool"] == "meeko" for t in res.report["tools"])


def test_meeko_takes_over_when_the_table_lacks_rules(tmp_path, monkeypatch):
    t = dict(recwriter._table())
    t["types"] = {k: v for k, v in t["types"].items() if k != "PHE"}
    monkeypatch.setattr(recwriter, "_table", lambda: t)
    p = tmp_path / "p.pdb"
    peptide_structure().write_pdb(str(p))
    res = prepare_receptor(p, tmp_path / "o", ReceptorPrepOptions())
    assert res.report["pdbqt"]["writer"] == "meeko" and res.report["status"] == "ok"
    assert any("no rules for" in w and "Meeko wrote the PDBQT" in w for w in res.report["warnings"])
    assert any(l[17:20] == "PHE" for l in open(res.pdbqt))


def test_table_only_reports_the_residues_it_cannot_type(tmp_path, monkeypatch):
    t = dict(recwriter._table())
    t["types"] = {k: v for k, v in t["types"].items() if k != "PHE"}
    monkeypatch.setattr(recwriter, "_table", lambda: t)
    p = tmp_path / "p.pdb"
    peptide_structure().write_pdb(str(p))
    res = prepare_receptor(p, tmp_path / "o", ReceptorPrepOptions(pdbqt_writer="table"))
    assert res.report["status"] == "degraded"
    assert any("PHE" in d["residue"] for d in res.report["dropped_residues"])


def test_compare_pdbqt_flags_a_type_disagreement_and_ignores_renames(tmp_path):
    from dockchem.recprep import compare_pdbqt

    def line(name, ad, res="ALA", num=1):
        return (f"ATOM  {1:5d} {name:<4s} {res:>3s} A{num:4d}       0.000   0.000   0.000  1.00  0.00    +0.000 {ad:<2s}\n")
    a, b, c = tmp_path / "a", tmp_path / "b", tmp_path / "c"
    a.write_text(line("N", "N") + line("H", "HD"))
    b.write_text(line("N", "N") + line("H1", "HD"))               # Meeko renames the N-terminal hydrogen
    c.write_text(line("N", "NA") + line("H", "HD"))
    assert compare_pdbqt(a, b)["result"] == "agree"
    r = compare_pdbqt(a, c)
    assert r["result"] == "disagree" and "A:ALA1" in r["examples"][0]
    d = tmp_path / "d"
    d.write_text(line("N", "N") + line("H", "HD") + line("OXT", "OA"))          # table wrote a terminal atom Meeko omitted
    e = compare_pdbqt(d, a)
    assert e["result"] == "agree" and e["n_residues_table_has_extra_atoms"] == 1
    assert compare_pdbqt(a, d)["result"] == "disagree"                      # Meeko wrote something the table lacks


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


# ---- kept hetero groups (cofactors, sugars) ---------------------------------------------------------------

from pathlib import Path as _P

DATA = _P(__file__).parent / "data"


def test_cofactor_is_typed_through_the_dictionary_route():
    """Real NADP block from 1T40: bond orders from the wwPDB dictionary, hydrogens added, Meeko atom types."""
    import collections
    atoms, warns = recwriter.type_hetero_group((DATA / "NAP_1T40.pdb").read_text(), "NAP", "A", 318)
    types = collections.Counter(a[2] for a in atoms)
    assert types["P"] == 3 and types["OA"] >= 10 and types["HD"] >= 5 and types["A"] >= 6
    assert all(t in {"C", "A", "N", "NA", "OA", "P", "S", "SA", "HD", "F", "Cl", "Br", "I"} for t in types)
    assert any("dictionary protonation state" in w for w in warns)
    assert len({a[0] for a in atoms}) == len(atoms)                       # unique atom names


def test_metal_cofactor_is_refused_with_the_reason():
    from dockchem.errors import HelperError
    with pytest.raises(HelperError, match="metal"):
        recwriter.type_hetero_group((DATA / "HEM_1MMV.pdb").read_text(), "HEM", "A", 1750)


def test_settle_pending_types_what_it_can_and_reports_the_rest(tmp_path):
    from dockchem.recprep import _Report, _settle_pending_het
    rep = _Report()
    rep.pending_het = {("A", "318", "NAP"): (DATA / "NAP_1T40.pdb").read_text().splitlines(),
                       ("A", "1750", "HEM"): (DATA / "HEM_1MMV.pdb").read_text().splitlines()}
    out = tmp_path / "r.pdbqt"
    out.write_text("")
    _settle_pending_het(rep, out, 0)
    text = out.read_text()
    assert " NAP A 318" in text and "HEM" not in text
    assert [d["residue"] for d in rep.dropped] == ["A:HEM1750"] and "metal" in rep.dropped[0]["reason"]
    assert any(c["stage"] == "hetero_groups" and "NAP" in c["residue"] for c in rep.changes)


def test_settle_pending_with_meeko_writer_reports_everything_as_left_out():
    from dockchem.recprep import _Report, _settle_pending_het
    rep = _Report()
    rep.pending_het = {("A", "318", "NAP"): (DATA / "NAP_1T40.pdb").read_text().splitlines()}
    _settle_pending_het(rep, None)
    assert rep.dropped and rep.dropped[0]["reason"] == "no bundled Meeko template"
