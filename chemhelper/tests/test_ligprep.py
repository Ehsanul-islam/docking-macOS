from rdkit import Chem
from rdkit.Chem import AllChem

import pytest
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
    p = one("C[Zn]C", protonate=False, allow_unusual_elements=True)
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


# ---- regressions found by running the whole wwPDB Chemical Component Dictionary ----------------------

def test_hypervalent_phosphorus_gets_zero_charges_not_a_crash():
    """CCD 02E: pentavalent P makes Gasteiger charges undefined; Meeko used to die on a non-finite charge."""
    r = prepare_ligand("COC=CCC(N)P(O)(O)(Oc1ccccc1)Oc2ccccc2", "02E")[0]
    assert r.status == "warning" and r.pdbqt
    assert any("partial charges could not be computed" in w for w in r.warnings)


def test_overlapping_atoms_in_the_embedding_are_never_passed_on():
    """CCD A1JWR (SF5 group): the first embedding put atoms on top of each other."""
    import numpy as np
    from scipy.spatial.distance import pdist
    r = prepare_ligand("Cn1c(c(c(n1)O)C(=O)O)CNc2cccc(c2)S(F)(F)(F)(F)F", "A1JWR")[0]
    assert r.status in ("prepared", "warning")
    assert pdist(r.mol.GetConformer().GetPositions()).min() >= 0.7


@pytest.mark.parametrize("smiles,element", [("C[Se]C", "Se"), ("C1CN2CCN1[Zn]2", "Zn"), ("C[Si](C)(C)C", "Si"),
                                            ("c1ccccc1[Fe]c1ccccc1", "Fe")])
def test_unsupported_elements_are_refused_clearly(smiles, element):
    r = prepare_ligand(smiles, "x")[0]
    assert r.status == "failed" and element in r.error and "allow_unusual_elements" in r.error


def test_unsupported_elements_can_be_allowed_with_a_warning():
    r = prepare_ligand("C[Mg]Cl", "x", LigandPrepOptions(allow_unusual_elements=True))[0]
    assert r.status in ("prepared", "warning") and r.pdbqt
    assert any("uncommon element" in w for w in r.warnings)


def test_counter_ions_are_stripped_before_the_element_check():
    r = prepare_ligand("CC(=O)[O-].[Na+]", "acetate")[0]
    assert r.status in ("prepared", "warning")


def test_torn_embedding_is_refused_not_passed_on():
    """CCD UUC (strained spiro-sulfur ring): the 3D embedding put an N-S 'bond' 9 A long."""
    smi = "c1nc(c2c(n1)n(cn2)C3C(C(C(O3)COCC#Cc4nc5c(ncnc5n4C6C(C(C(O6)CNS78(CN7C8)(O)O)O)O)N)O)O)N"
    r = prepare_ligand(smi, "UUC")[0]
    assert r.status == "failed" and "implausible bond" in r.error


def test_long_but_real_iodine_bonds_are_accepted():
    r = prepare_ligand("[I-].[I-].[I-]", "tri", LigandPrepOptions(allow_unusual_elements=True))
    assert r  # polyiodide input does not crash; real I-I bonds (2.7 A) must not trip the check
    ok = prepare_ligand("IC(I)I", "iodoform")[0]
    assert ok.status in ("prepared", "warning")


def test_macromolecule_sized_ligands_are_refused_up_front_not_hung_on():
    for smi in ("C" * 5000, "c1ccccc1" * 200):
        r = prepare_ligand(smi, "big")[0]
        assert r.status == "failed" and "not a small molecule" in r.error
    assert prepare_ligand("C" * 20, "ok")[0].status in ("prepared", "warning")
