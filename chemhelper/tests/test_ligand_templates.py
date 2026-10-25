"""Bond orders for ligands cut from structure files. Uses the real PDB entries in pdb_structures_50/."""
import json
import subprocess
import sys
from pathlib import Path

import pytest
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

from dockchem import formats, ligand_templates as lt
from dockchem.errors import InvalidInput, NeedsTemplate
from dockchem.ligprep import prepare_ligand

PDBS = Path(__file__).resolve().parents[2] / "pdb_structures_50"
pytestmark = pytest.mark.skipif(not PDBS.is_dir(), reason="pdb_structures_50 not present")


def pdb(pid):
    return PDBS / f"{pid}.pdb"


def formula(mol):
    return rdMolDescriptors.CalcMolFormula(mol)


# ---- the audit finding: imatinib came out saturated, biotin's acid became a gem-diol -----------------
def test_imatinib_from_1IEP_has_correct_chemistry():
    r = lt.ligand_from_structure(pdb("1IEP"), "STI")
    assert formula(r.mol) == "C29H31N7O" and r.source == "dictionary:STI"
    flat = Chem.RemoveHs(r.mol)
    assert sum(a.GetIsAromatic() for a in flat.GetAtoms()) == 24          # four six-membered aromatic rings, not saturated
    assert flat.HasSubstructMatch(Chem.MolFromSmarts("c1ccccc1NC(=O)c1ccccc1"))
    assert r.audit == []


def test_biotin_from_1STP_keeps_its_acid_and_urea():
    r = lt.ligand_from_structure(pdb("1STP"), "BTN")
    assert formula(r.mol) == "C10H16N2O3S"
    assert r.mol.HasSubstructMatch(Chem.MolFromSmarts("C(=O)[OH]"))
    assert r.mol.HasSubstructMatch(Chem.MolFromSmarts("NC(=O)N"))
    assert len(Chem.FindMolChiralCenters(r.mol, useLegacyImplementation=False)) == 3     # from the crystal


def test_no_phantom_phosphorus_stereo_on_atp():
    r = lt.ligand_from_structure(pdb("1ATP"), "ATP")
    assert all(a.GetChiralTag() == Chem.ChiralType.CHI_UNSPECIFIED for a in r.mol.GetAtoms() if a.GetSymbol() == "P")


# ---- the safety net --------------------------------------------------------------------------------
def test_ligand_read_from_pdb_is_untrusted_and_refused_everywhere(tmp_path):
    lig = tmp_path / "ligand.pdb"
    mol, _, _ = lt.residue_to_mol(pdb("1IEP"), "STI")
    Chem.MolToPDBFile(mol, str(lig))
    mols, warns = formats.read_molecules(lig)
    assert lt.is_untrusted(mols[0]) and any("heavy-atom skeleton" in w for w in warns)
    r = prepare_ligand(mols[0], "x")[0]
    assert r.status == "failed" and "bond orders are unknown" in r.error
    for ext in (".sdf", ".mol2", ".smi"):
        with pytest.raises(NeedsTemplate):
            formats.write_molecules(mols, tmp_path / f"o{ext}")
    with pytest.raises(NeedsTemplate):
        formats.convert(lig, tmp_path / "o.sdf")


def test_unknown_ligand_asks_for_smiles_instead_of_guessing():
    with pytest.raises(NeedsTemplate) as e:
        lt.ligand_from_structure(pdb("4G5J"), "0WN")
    assert e.value.needs == "smiles" and e.value.resname == "0WN" and "not in the bundled" in e.value.message


def test_wrong_template_is_rejected_not_forced():
    with pytest.raises(NeedsTemplate):                            # different composition entirely
        lt.ligand_from_structure(pdb("3PTB"), "BEN", smiles="CCO")
    with pytest.raises(NeedsTemplate):                            # same composition, different graph
        lt.ligand_from_structure(pdb("3PTB"), "BEN", smiles="CC(C)(C)C#N")


def test_user_smiles_overrides_dictionary_and_allows_charged_form():
    r = lt.ligand_from_structure(pdb("3PTB"), "BEN", smiles="NC(=[NH2+])c1ccccc1")
    assert r.source == "user_smiles" and Chem.GetFormalCharge(r.mol) == 1


def test_unparseable_smiles_is_an_input_error():
    with pytest.raises(InvalidInput):
        lt.ligand_from_structure(pdb("3PTB"), "BEN", smiles="not a smiles")


def test_missing_residue_is_an_input_error():
    with pytest.raises(InvalidInput):
        lt.ligand_from_structure(pdb("3PTB"), "ZZZ")


def test_net_charge_works_only_when_hydrogens_are_present():
    r = lt.ligand_from_structure(pdb("7RSA"), "TBU", charge=0)        # neutron-style entry with hydrogens
    assert formula(r.mol) == "C4H10O" and r.source == "perceived_with_charge"
    with pytest.raises(NeedsTemplate):
        lt.ligand_from_structure(pdb("7RSA"), "TBU", charge=1)        # wrong charge must not be accepted
    with pytest.raises(NeedsTemplate) as e:
        lt.ligand_from_structure(pdb("4G5J"), "0WN", charge=0)        # no hydrogens: charge alone is not enough
    assert "no hydrogens" in e.value.message


# ---- dictionary integrity --------------------------------------------------------------------------
def test_dictionary_is_nonempty_and_every_entry_names_its_validation():
    d = lt.load_dictionary()
    assert len(d) >= 10
    for res, e in d.items():
        assert e["validated_on"] and Chem.MolFromSmiles(e["smiles"]) is not None, res


@pytest.mark.parametrize("res", sorted(lt.load_dictionary()))
def test_every_dictionary_entry_still_matches_its_crystal_ligand(res):
    e = lt.load_dictionary()[res]
    for pid in e["validated_on"]:
        r = lt.ligand_from_structure(pdb(pid), res)
        assert r.source == f"dictionary:{res}" and r.audit == [], (res, pid)


def test_dictionary_ships_inside_the_package():
    assert lt.DICT_PATH.is_file() and lt.DICT_PATH.parent.name == "data"
    assert json.loads(lt.DICT_PATH.read_text())["version"]


# ---- end to end ------------------------------------------------------------------------------------
def test_resolved_ligand_prepares_for_docking_and_keeps_crystal_coordinates():
    r = lt.ligand_from_structure(pdb("1IEP"), "STI")
    ref = r.mol.GetConformer().GetPositions().copy()
    p = prepare_ligand(r.mol, "STI")[0]
    assert p.status in ("prepared", "warning") and "ROOT" in p.pdbqt, p.error
    assert abs(p.mol.GetConformer().GetPositions() - ref).max() < 1e-4        # crystal pose preserved


def test_receptor_prep_reports_bond_order_status_for_extracted_ligands(tmp_path):
    from dockchem.recprep import ReceptorPrepOptions, prepare_receptor
    ok = prepare_receptor(pdb("1STP"), tmp_path / "a", ReceptorPrepOptions(reference_ligand="BTN"))
    bo = ok.ligands[0]["bond_orders"]
    assert bo["status"] == "resolved" and Path(bo["sdf"]).exists()
    mol = next(Chem.SDMolSupplier(bo["sdf"], removeHs=False))
    assert formula(mol) == "C10H16N2O3S"
    bad = prepare_receptor(pdb("4G5J"), tmp_path / "b", ReceptorPrepOptions(reference_ligand="0WN"))
    assert bad.ligands[0]["bond_orders"]["status"] == "needs_template"


def test_rpc_from_structure_and_typed_error(tmp_path):
    def call(*reqs):
        p = subprocess.run([sys.executable, "-m", "dockchem.server"], capture_output=True, text=True, timeout=120,
                           input="\n".join(json.dumps({"jsonrpc": "2.0", "id": i, "method": m, "params": pr})
                                           for i, (m, pr) in enumerate(reqs, 1)) + "\n")
        return {m["id"]: m for m in map(json.loads, p.stdout.splitlines())}
    out = call(("ligand.from_structure", {"path": str(pdb("1STP")), "resname": "BTN", "outdir": str(tmp_path)}),
               ("ligand.from_structure", {"path": str(pdb("4G5J")), "resname": "0WN"}))
    assert out[1]["result"]["result"]["heavy_atoms"] == 16 and (tmp_path / "BTN.sdf").exists()
    err = out[2]["error"]["data"]
    assert err["type"] == "NeedsTemplate" and err["needs"] == "smiles" and err["resname"] == "0WN"
