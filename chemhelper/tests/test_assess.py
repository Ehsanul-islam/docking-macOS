"""structure.assess: every verdict, and the gate in prepare_receptor."""
import copy

import gemmi
import pytest

from dockchem import methods
from dockchem.assess import assess, assess_structure
from dockchem.errors import InvalidInput, UnusableStructure
from dockchem.recprep import ReceptorPrepOptions, prepare_receptor

from fixtures import peptide_structure

SEQ30 = "GSHDEKRCNQFYTALGSHDEKRCNQFYTAL"


@pytest.fixture(scope="module")
def base():
    return peptide_structure(SEQ30)


def codes(a, sev=None):
    return {f["code"] for f in a["findings"] if sev is None or f["severity"] == sev}


def write(st, tmp_path, name="s.pdb"):
    p = tmp_path / name
    st.write_pdb(str(p))
    return p


def test_clean_structure_is_ready(base, tmp_path):
    a = assess(write(base, tmp_path))
    assert a["verdict"] == "ready" and a["summary"] == "no problems found"
    assert a["stats"]["protein_residues"] == 30 and a["stats"]["distorted_residues"] == 0


def test_ca_only_model_is_unusable(base, tmp_path):
    st = copy.deepcopy(base)
    for r in st[0][0]:
        for i in reversed(range(len(r))):
            if r[i].name != "CA":
                del r[i]
    a = assess(write(st, tmp_path))
    assert a["verdict"] == "unusable" and "ca_only" in codes(a, "block")
    assert "Calpha-only" in a["summary"]


def test_backbone_only_model_is_unusable(base, tmp_path):
    st = copy.deepcopy(base)
    for r in st[0][0]:
        for i in reversed(range(len(r))):
            if r[i].name not in ("N", "CA", "C", "O"):
                del r[i]
    assert "backbone_only" in codes(assess(write(st, tmp_path)), "block")


def test_nucleic_acid_only_is_unusable(tmp_path):
    p = tmp_path / "dna.pdb"
    lines = []
    for i, base_ in enumerate(("DA", "DT", "DG", "DC"), 1):
        for j, (n, el) in enumerate((("P", "P"), ("O5'", "O"), ("C5'", "C"), ("C4'", "C"))):
            lines.append(f"ATOM  {i * 10 + j:5d} {n:<4s} {base_:>3s} A{i:4d}    {i * 3.4:8.3f}{j * 1.4:8.3f}   0.000  1.00  0.00          {el:>2s}")
    p.write_text("\n".join(lines) + "\nEND\n")
    a = assess(p)
    assert a["verdict"] == "unusable" and "nucleic_acid_only" in codes(a)


@pytest.mark.parametrize("content", ["", "this is not a structure\n", "ATOM garbage\nEND\n"])
def test_garbage_files_get_a_verdict_not_a_crash(tmp_path, content):
    p = tmp_path / "x.pdb"
    p.write_text(content)
    assert assess(p)["verdict"] == "unusable"


def test_missing_file_is_unusable(tmp_path):
    assert assess(tmp_path / "nope.pdb")["verdict"] == "unusable"


def test_missing_side_chain_atoms_are_fixable_and_counted(base, tmp_path):
    st = copy.deepcopy(base)
    r = [x for x in st[0][0] if x.name == "LYS"][0]
    for i in reversed(range(len(r))):
        if r[i].name in ("CE", "NZ"):
            del r[i]
    a = assess(write(st, tmp_path))
    assert a["verdict"] == "fixable" and "missing_atoms" in codes(a, "fix")
    assert a["stats"]["residues_with_missing_atoms"] == 1


def test_one_distorted_residue_is_fixable_many_are_unusable(base, tmp_path):
    one = copy.deepcopy(base)
    r = one[0][0][5]
    r["CA"][0].pos = gemmi.Position(r["N"][0].pos.x + 0.5, r["N"][0].pos.y, r["N"][0].pos.z)
    a = assess(write(one, tmp_path, "one.pdb"))
    assert a["verdict"] == "fixable" and "distorted_residues" in codes(a, "fix")
    many = copy.deepcopy(base)
    for r in list(many[0][0])[:6]:
        r["CA"][0].pos = gemmi.Position(r["N"][0].pos.x + 0.5, r["N"][0].pos.y, r["N"][0].pos.z)
    b = assess(write(many, tmp_path, "many.pdb"))
    assert b["verdict"] == "unusable" and "distorted_coordinates" in codes(b, "block")


def test_chain_break_is_a_note_only(base, tmp_path):
    st = copy.deepcopy(base)
    for r in list(st[0][0])[15:]:
        for a in r:
            a.pos = gemmi.Position(a.pos.x + 30, a.pos.y, a.pos.z)
    a = assess(write(st, tmp_path))
    assert a["verdict"] == "ready" and a["stats"]["chain_breaks"] == 1 and "chain_breaks" in codes(a, "note")


def test_alternate_conformations_are_fixable(base, tmp_path):
    st = copy.deepcopy(base)
    for a in st[0][0][3]:
        a.altloc = "A"
    assert "alternate_locations" in codes(assess(write(st, tmp_path)), "fix")


def test_unknown_residues(base, tmp_path):
    some = copy.deepcopy(base)
    some[0][0][4].name = "UNK"
    assert "unknown_residues" in codes(assess(write(some, tmp_path, "a.pdb")), "fix")
    most = copy.deepcopy(base)
    for r in list(most[0][0])[:8]:
        r.name = "UNK"
    assert "unknown_residues" in codes(assess(write(most, tmp_path, "b.pdb")), "block")


def test_multi_model_noted_and_bad_model_number(base, tmp_path):
    st = copy.deepcopy(base)
    second = copy.deepcopy(st[0])
    second.num = 2
    st.add_model(second)
    p = write(st, tmp_path)
    assert "multiple_models" in codes(assess(p))
    with pytest.raises(InvalidInput):
        assess(p, model=5)


def test_only_peptide_note_for_short_chains(tmp_path):
    a = assess(write(peptide_structure(), tmp_path))
    assert "short_chain" in codes(a, "note") and a["verdict"] == "ready"


def test_prepare_refuses_unusable_with_reasons_and_allows_override(base, tmp_path):
    bad = copy.deepcopy(base)
    for r in list(bad[0][0])[:6]:
        r["CA"][0].pos = gemmi.Position(r["N"][0].pos.x + 0.5, r["N"][0].pos.y, r["N"][0].pos.z)
    p = write(bad, tmp_path, "bad.pdb")
    with pytest.raises(UnusableStructure, match="distorted coordinates") as ei:
        prepare_receptor(p, tmp_path / "o", ReceptorPrepOptions())
    rpc = ei.value.to_rpc()
    assert rpc["data"]["assessment"]["verdict"] == "unusable"
    res = prepare_receptor(p, tmp_path / "o2", ReceptorPrepOptions(allow_unusable=True))
    assert any("prepared on request" in w for w in res.report["warnings"])
    assert res.report["assessment"]["verdict"] == "unusable"


def test_assessment_is_in_the_prep_report_and_edits_can_fix_a_bad_chain(base, tmp_path):
    ok = prepare_receptor(write(base, tmp_path), tmp_path / "o", ReceptorPrepOptions())
    assert ok.report["assessment"]["verdict"] == "ready"
    st = copy.deepcopy(base)
    ch = gemmi.Chain("B")
    for r in copy.deepcopy(list(st[0][0])):
        for a in r:
            a.pos = gemmi.Position(a.pos.x + 40, a.pos.y, a.pos.z)
            a.pos = a.pos
        ch.add_residue(r)
    for r in list(ch)[:6]:
        r["CA"][0].pos = gemmi.Position(r["N"][0].pos.x + 0.5, r["N"][0].pos.y, r["N"][0].pos.z)
    st[0].add_chain(ch)
    p = write(st, tmp_path, "two.pdb")
    assert assess(p)["verdict"] in ("fixable", "unusable")
    res = prepare_receptor(p, tmp_path / "o2", ReceptorPrepOptions(edits=["delete chain B"]))
    assert res.report["assessment"]["verdict"] == "ready"


def test_rpc_method(base, tmp_path):
    out = methods.REGISTRY["structure.assess"]({"path": str(write(base, tmp_path))}, None) \
        if hasattr(methods, "REGISTRY") else methods.structure_assess({"path": str(write(base, tmp_path))}, None)
    assert out["result"]["verdict"] == "ready" and out["provenance"]["tool"] == "dockchem.assess"


def test_short_x_ray_hydroxyl_hydrogens_are_not_distortion(base, tmp_path):
    """Refined entries (5SAK, 6XCT, ...) put OG-HG at 0.84 A; 49 of 300 benchmark entries were wrongly refused."""
    st = copy.deepcopy(base)
    n = 0
    for r in st[0][0]:
        if r.name in ("SER", "THR", "TYR"):
            o = next(a for a in r if a.name in ("OG", "OG1", "OH"))
            h = gemmi.Atom()
            h.name, h.element, h.pos = "H" + o.name[1:], gemmi.Element("H"), gemmi.Position(o.pos.x + 0.84, o.pos.y, o.pos.z)
            r.add_atom(h)
            n += 1
    assert n >= 3
    a = assess(write(st, tmp_path))
    assert a["stats"]["distorted_residues"] == 0 and "distorted_coordinates" not in codes(a)
