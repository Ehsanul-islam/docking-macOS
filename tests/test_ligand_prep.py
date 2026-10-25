import io
import json
from pathlib import Path

import pytest
from rdkit import Chem

from dockmac import ligand_prep as L
from dockmac import structure as S

PDB = Path(__file__).parent.parent / "pdb_structures_50"
ERLOTINIB = "COCCOc1cc2ncnc(Nc3cccc(C#C)c3)c2cc1OCCOC"


def test_smiles_gets_hydrogens_and_3d():
    mol = L.from_smiles("CCO")
    assert mol.GetNumAtoms() == 9
    assert mol.GetConformer().Is3D()


def test_invalid_smiles():
    with pytest.raises(ValueError):
        L.from_smiles("not a smiles")


def test_cocrystal_ligand_with_template():
    st = S.load(PDB / "1M17.pdb")
    lig = S.extract_ligand(st, "AQ4")
    mol = L.from_structure(lig, smiles_template=ERLOTINIB)
    assert Chem.MolToSmiles(Chem.RemoveHs(mol)) == Chem.MolToSmiles(Chem.MolFromSmiles(ERLOTINIB))
    assert sum(1 for a in mol.GetAtoms() if a.GetSymbol() == "H") == 23


def test_template_mismatch_raises():
    lig = S.extract_ligand(S.load(PDB / "1M17.pdb"), "AQ4")
    with pytest.raises(ValueError):
        L.from_structure(lig, smiles_template="CCO")


def test_fetch_template_smiles(monkeypatch):
    class Resp:
        def __enter__(self):
            return io.BytesIO(json.dumps({"rcsb_chem_comp_descriptor": {"smiles": ERLOTINIB}}).encode())

        def __exit__(self, *a):
            pass

    monkeypatch.setattr(L.urllib.request, "urlopen", lambda *a, **k: Resp())
    assert L.fetch_template_smiles("aq4") == ERLOTINIB


def test_pdbqt_output():
    text = L.to_pdbqt(L.from_smiles(ERLOTINIB))
    assert text.startswith("ROOT") or "ROOT" in text
    assert "TORSDOF" in text
