from pathlib import Path

from dockmac import protein_prep as P
from dockmac import receptor_pdbqt as R
from dockmac import structure as S

PDB = Path(__file__).parent.parent / "pdb_structures_50"


def _prepared(name):
    st = S.clean(S.load(PDB / f"{name}.pdb"), remove_ligands=True)
    return P.prepare_protein(st, ph=7.4)[0]


def test_lysozyme_net_charge_and_types():
    text, rep = R.to_pdbqt(_prepared("1LYZ"))
    atoms = [l for l in text.splitlines() if l.startswith("ATOM")]
    assert rep.charges_from_forcefield and rep.n_atoms == len(atoms)
    assert round(sum(float(l[70:76]) for l in atoms)) == 8
    types = {l[77:79].strip() for l in atoms}
    assert {"C", "A", "N", "OA", "HD"} <= types
    assert "H" not in types  # non-polar hydrogens merged


def test_every_pdb_gets_a_receptor():
    for name in ["1UBQ", "1M17", "2RH1"]:
        text, _ = R.to_pdbqt(_prepared(name))
        assert text.endswith("END\n")
