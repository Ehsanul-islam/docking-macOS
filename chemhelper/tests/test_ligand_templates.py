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
    assert formula(r.mol) == "C29H31N7O" and r.source == "ccd:STI"
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


def test_ligand_absent_from_every_dictionary_asks_for_smiles(tmp_path):
    f = tmp_path / "x.pdb"
    assert lt.ccd_entry("X0Z") is None                                      # a name the CCD does not have
    f.write_text(pdb("1STP").read_text().replace(" BTN ", " X0Z "))
    with pytest.raises(NeedsTemplate) as e:
        lt.ligand_from_structure(f, "X0Z")
    assert e.value.needs == "smiles" and e.value.resname == "X0Z" and "not in the bundled" in e.value.message
    ok = lt.ligand_from_structure(f, "X0Z", smiles="OC(=O)CCCCC1SCC2NC(=O)NC12")      # the user's way out
    assert ok.source == "user_smiles" and formula(ok.mol) == "C10H16N2O3S"


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
    mol, _, _ = lt.residue_to_mol(pdb("7RSA"), "TBU")                 # neutron-style entry with hydrogens
    r = lt._perceive_with_charge(mol, 0, "TBU")
    assert formula(r) == "C4H10O" and r.GetProp(lt.SOURCE) == "perceived_with_charge"
    with pytest.raises(NeedsTemplate):
        lt._perceive_with_charge(mol, 1, "TBU")                       # wrong charge must not be accepted
    mol, _, _ = lt.residue_to_mol(pdb("4G5J"), "0WN")
    with pytest.raises(NeedsTemplate) as e:
        lt._perceive_with_charge(mol, 0, "0WN")                        # no hydrogens: charge alone is not enough
    assert "no hydrogens" in e.value.message


# ---- the Chemical Component Dictionary bundle -------------------------------------------------------
# Ligands that must keep resolving, with the structure that proves it (all checked against the CCD).
REFERENCE = {"STI": "1IEP", "AQ4": "1M17", "BTN": "1STP", "BEN": "3PTB", "CAU": "2RH1", "ATP": "1ATP",
             "GTP": "1JFF", "GDP": "1JFF", "NAD": "4MDH", "NDP": "8CAT", "AP5": "1AKE", "CLR": "2RH1",
             "BU1": "2RH1", "0WN": "4G5J", "0WM": "4G5J", "F09": "1K4C", "P0G": "3SN6", "TXL": "1TUB",
             "TA1": "1JFF"}


def test_ccd_bundle_is_present_complete_and_dated():
    info = lt.ccd_info()
    assert info["count"] > 50000 and "wwPDB" in info["source"] and info["last_modified"]
    assert lt.ccd_entry("ATP")["f"] == "C10 H16 N5 O13 P3" and lt.ccd_entry("atp") == lt.ccd_entry("ATP")
    assert lt.CCD_PATH.stat().st_size < 10_000_000                     # compact enough to ship


@pytest.mark.parametrize("res,pid", sorted(REFERENCE.items()))
def test_reference_ligands_resolve_from_the_ccd_with_a_clean_audit(res, pid):
    r = lt.ligand_from_structure(pdb(pid), res)
    assert r.source == f"ccd:{res}" and r.audit == [], (res, pid, r.audit)
    assert any("Chemical Component Dictionary" in w for w in r.warnings)


def test_metal_complexes_are_refused_with_a_clear_reason():
    with pytest.raises(NeedsTemplate) as e:
        lt.ligand_from_structure(pdb("1MBN"), "HEM")
    assert "metal complex" in e.value.message and e.value.needs == "none"


def test_truncated_or_covalent_ligands_are_refused_not_capped():
    for pid, res in (("2RH1", "PLM"), ("1HZH", "NAG"), ("1K4C", "DGA")):      # atoms missing from the file
        with pytest.raises(NeedsTemplate) as e:
            lt.ligand_from_structure(pdb(pid), res)
        assert "does not fit the atoms in this file" in e.value.message, res


def test_ccd_entry_agrees_with_independent_hand_written_smiles():
    """Cross-check done once by hand for 13 ligands (identical results); kept as a regression for three."""
    for res, smi in (("STI", "Cc1ccc(NC(=O)c2ccc(CN3CCN(C)CC3)cc2)cc1Nc1nccc(n1)-c1cccnc1"),
                     ("BTN", "OC(=O)CCCCC1SCC2NC(=O)NC12"), ("CAU", "CC(C)NCC(O)COc1cccc2[nH]c3ccccc3c12")):
        pid = REFERENCE[res]
        a = lt.ligand_from_structure(pdb(pid), res)
        b = lt.ligand_from_structure(pdb(pid), res, smiles=smi)
        assert Chem.MolToSmiles(Chem.RemoveHs(a.mol)) == Chem.MolToSmiles(Chem.RemoveHs(b.mol))


def test_distiller_parses_a_real_ccd_block():
    import importlib.util
    spec = importlib.util.spec_from_file_location(
        "build_ccd", Path(__file__).resolve().parents[2] / "tools" / "ligand_dictionary" / "build_ccd.py")
    mod = importlib.util.module_from_spec(spec); spec.loader.exec_module(mod)
    cid, e = mod.distil((Path(__file__).parent / "data" / "BTN.cif").read_text())
    assert cid == "BTN" and e["f"] == "C10 H16 N2 O3 S" and e["n"] == "BIOTIN"
    assert e["s"][0] == "C1C2C(C(S1)CCCCC(=O)O)NC(=O)N2"               # flat OpenEye SMILES preferred


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
    bad = prepare_receptor(pdb("1MBN"), tmp_path / "b", ReceptorPrepOptions(reference_ligand="HEM"))
    assert bad.ligands[0]["bond_orders"]["status"] == "needs_template"           # heme: metal complex
    assert "metal complex" in bad.ligands[0]["bond_orders"]["message"]


def test_rpc_from_structure_and_typed_error(tmp_path):
    def call(*reqs):
        p = subprocess.run([sys.executable, "-m", "dockchem.server"], capture_output=True, text=True, timeout=120,
                           input="\n".join(json.dumps({"jsonrpc": "2.0", "id": i, "method": m, "params": pr})
                                           for i, (m, pr) in enumerate(reqs, 1)) + "\n")
        return {m["id"]: m for m in map(json.loads, p.stdout.splitlines())}
    out = call(("ligand.from_structure", {"path": str(pdb("1STP")), "resname": "BTN", "outdir": str(tmp_path)}),
               ("ligand.from_structure", {"path": str(pdb("1MBN")), "resname": "HEM"}))
    assert out[1]["result"]["result"]["heavy_atoms"] == 16 and (tmp_path / "BTN.sdf").exists()
    err = out[2]["error"]["data"]
    assert err["type"] == "NeedsTemplate" and err["needs"] == "none" and err["resname"] == "HEM"
