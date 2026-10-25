from pathlib import Path

from dockmac import protein_prep as P
from dockmac import structure as S

PDB = Path(__file__).parent.parent / "pdb_structures_50"


def _heavy(st):
    return sum(1 for ch in st[0] for r in ch for a in r if not a.is_hydrogen())


def _h(st):
    return sum(1 for ch in st[0] for r in ch for a in r if a.is_hydrogen())


def test_prepare_adds_hydrogens_and_drops_hetero():
    st = S.clean(S.load(PDB / "1UBQ.pdb"))
    out, rep = P.prepare_protein(st)
    assert rep.hydrogens_added > 600 and _h(out) == rep.hydrogens_added
    assert not any(h for c in S.summarize(out) for h in c.hetero)


def test_prepare_does_not_modify_input():
    st = S.clean(S.load(PDB / "1UBQ.pdb"))
    n = st[0].count_atom_sites()
    P.prepare_protein(st)
    assert st[0].count_atom_sites() == n


def test_no_hydrogens_option():
    st = S.clean(S.load(PDB / "1UBQ.pdb"))
    out, rep = P.prepare_protein(st, add_hydrogens=False)
    assert _h(out) == 0 and rep.hydrogens_added == 0


def test_receptor_with_ligand_is_prepared():
    st = S.clean(S.load(PDB / "1M17.pdb"), remove_ligands=True)
    out, rep = P.prepare_protein(st, ph=7.4)
    assert _heavy(out) > 2000 and rep.hydrogens_added > 2000
