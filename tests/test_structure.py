from pathlib import Path

import pytest

from dockmac import structure as S

PDB = Path(__file__).parent.parent / "pdb_structures_50"


def test_summarize_1m17_finds_erlotinib():
    st = S.load(PDB / "1M17.pdb")
    info = S.summarize(st)
    ligands = {h.resname for c in info for h in c.hetero if h.kind == "ligand"}
    assert ligands == {"AQ4"}  # erlotinib
    assert sum(c.n_polymer_residues for c in info) > 200


def test_clean_removes_water_and_ions_keeps_ligand():
    st = S.load(PDB / "1M17.pdb")
    out = S.clean(st)
    kinds = {h.kind for c in S.summarize(out) for h in c.hetero}
    assert "water" not in kinds and "ion" not in kinds
    assert "ligand" in kinds


def test_clean_does_not_modify_input():
    st = S.load(PDB / "1M17.pdb")
    n = st[0].count_atom_sites()
    S.clean(st, remove_ligands=True)
    assert st[0].count_atom_sites() == n


def test_remove_ligands_and_keep_chain():
    st = S.load(PDB / "2HHB.pdb")  # hemoglobin, 4 chains, HEM
    out = S.clean(st, keep_chains=["A"], remove_ligands=True)
    info = S.summarize(out)
    assert [c.name for c in info] == ["A"]
    assert not any(h.kind == "ligand" for h in info[0].hetero)


def test_keep_ligand_by_name():
    st = S.load(PDB / "2HHB.pdb")
    out = S.clean(st, remove_ligands=True, keep_ligands=["HEM"])
    assert any(h.resname == "HEM" for c in S.summarize(out) for h in c.hetero)


def test_extract_ligand_and_missing():
    st = S.load(PDB / "2HHB.pdb")
    lig = S.extract_ligand(st, "HEM")
    assert lig[0].count_atom_sites() == 43
    with pytest.raises(ValueError):
        S.extract_ligand(st, "XXX")


def test_best_altloc_leaves_no_altlocs():
    for pdb in sorted(PDB.glob("*.pdb")):
        out = S.clean(S.load(pdb))
        assert not any(a.has_altloc() for ch in out[0] for r in ch for a in r), pdb.name
