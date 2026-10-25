from rdkit import Chem
from rdkit.Chem import AllChem

from dockchem.ligprep import (LigandPrepOptions, prepare_ligand, prepare_library, properties,
                              read_records)


def one(smi, **kw):
    r = prepare_ligand(smi, "x", LigandPrepOptions(**kw))
    assert len(r) == 1, [p.summary() for p in r]
    return r[0]


def test_basic_prep_produces_pdbqt_and_props():
    p = one("CC(=O)Nc1ccccc1")
    assert p.status == "prepared"
    assert "ROOT" in p.pdbqt and p.pdbqt.count("ATOM") >= 10
    assert p.properties["heavy_atoms"] == 10 and p.properties["lipinski_pass"]
    assert p.mol.GetNumConformers() == 1


def test_deterministic_with_seed():
    a, b = one("CCCCOc1ccccc1N", seed=7), one("CCCCOc1ccccc1N", seed=7)
    assert a.pdbqt == b.pdbqt


def test_salt_stripped_and_logged():
    p = one("CC(=O)Nc1ccc(O)cc1.[Na+].[Cl-]")
    assert any("largest of 3" in s for s in p.steps)
    assert p.smiles == one("CC(=O)Nc1ccc(O)cc1").smiles


def test_protonation_ph74():
    p = one("CC(=O)O")           # acetic acid -> carboxylate at 7.4
    assert "[O-]" in p.smiles
    q = one("CC(=O)O", protonate=False)
    assert "[O-]" not in q.smiles


def test_undefined_stereo_warns_or_enumerates():
    keep = one("CC(F)Cl", protonate=False)
    assert any("undefined stereocentre" in w for w in keep.warnings) and keep.status == "warning"
    isos = prepare_ligand("CC(F)Cl", "x", LigandPrepOptions(stereo="enumerate", protonate=False))
    assert len(isos) == 2 and {i.smiles for i in isos} == {"C[C@H](F)Cl", "C[C@@H](F)Cl"}


def test_defined_stereo_preserved():
    p = one("C[C@H](F)Cl", protonate=False)
    assert p.smiles == "C[C@H](F)Cl" and not p.warnings


def test_enumeration_cap():
    isos = prepare_ligand("CC(F)C(Cl)C(Br)C(I)C(O)C", "x",
                          LigandPrepOptions(stereo="enumerate", max_variants=4, protonate=False))
    assert len(isos) <= 4 and any("capped" in w for w in isos[0].warnings)


def test_input_3d_preserved():
    m = Chem.AddHs(Chem.MolFromSmiles("CC(=O)Nc1ccc(O)cc1"))
    AllChem.EmbedMolecule(m, randomSeed=3)
    ref = m.GetConformer().GetPositions()[:m.GetNumHeavyAtoms()].copy()
    p = prepare_ligand(m, "ref")[0]
    assert any("preserved input coordinates" in s for s in p.steps)
    out = p.mol.GetConformer().GetPositions()
    heavy = [a.GetIdx() for a in p.mol.GetAtoms() if a.GetAtomicNum() > 1]
    assert abs(out[heavy] - ref).max() < 1e-6


def test_bad_input_fails_cleanly():
    r = prepare_ligand("this is not smiles", "bad")[0]
    assert r.status == "failed" and r.error


def test_uncommon_element_warns():
    p = one("OB(O)c1ccccc1", protonate=False)
    assert any("uncommon element" in w for w in p.warnings)


def test_many_torsions_warns():
    p = one("C" * 20 + "O", protonate=False)
    assert any("torsions" in w for w in p.warnings)


def test_alerts_flagged_not_filtered():
    props = properties(Chem.AddHs(Chem.MolFromSmiles("O=C1C=CC(=O)C=C1")))   # quinone
    assert props["alerts"]


def test_library_mode(tmp_path):
    f = tmp_path / "lib.smi"
    f.write_text("CCO ethanol\nOCC ethanol_dup\nnot_smiles bad\nCC(=O)Nc1ccc(O)cc1 apap\n")
    lib = prepare_library(read_records(f), LigandPrepOptions(protonate=False))
    st = {r["name"]: r["status"] for r in lib["rows"]}
    assert st["ethanol"] in ("prepared", "warning") and st["ethanol_dup"] == "duplicate"
    assert st["bad"] == "failed" and st["apap"] == "prepared"
    assert lib["counts"]["failed"] == 1 and len(lib["prepared"]) == 2


def test_library_sdf_unreadable_record_listed(tmp_path):
    m = Chem.AddHs(Chem.MolFromSmiles("CCO")); AllChem.EmbedMolecule(m, randomSeed=1)
    w = Chem.SDWriter(str(tmp_path / "a.sdf")); m.SetProp("_Name", "eth"); w.write(m); w.close()
    recs = list(read_records(tmp_path / "a.sdf"))
    assert recs[0][0] == "eth"


def test_borderline_phenol_not_ionised_by_default():
    p = one("CC(=O)Nc1ccc(O)cc1")
    assert p.smiles == "CC(=O)Nc1ccc(O)cc1" and any("ambiguous" in w for w in p.warnings)


def test_ph_range_enumerates_states():
    r = prepare_ligand("CC(=O)Nc1ccc(O)cc1", "x", LigandPrepOptions(ph_range=2.0))
    assert len(r) == 2
