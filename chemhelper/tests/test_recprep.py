import copy

import pytest

from dockchem.errors import InvalidInput
from dockchem.recprep import ReceptorPrepOptions, prepare_receptor

from fixtures import decorated_structure, peptide_structure


@pytest.fixture(scope="module")
def decorated(tmp_path_factory):
    p = tmp_path_factory.mktemp("in") / "dec.pdb"
    decorated_structure().write_pdb(str(p))
    return p


@pytest.fixture(scope="module")
def prepared(decorated, tmp_path_factory):
    return prepare_receptor(decorated, tmp_path_factory.mktemp("out"),
                            ReceptorPrepOptions(reference_ligand="LIG"))


def changes(res, stage):
    return [c for c in res.report["changes"] if c["stage"] == stage]


def test_produces_pdbqt_with_hydrogens(prepared):
    txt = prepared.pdbqt.read_text()
    assert txt.count("ATOM") > 150 and " HD" in txt        # polar hydrogens kept, non-polar merged
    assert prepared.report["site_residues"]
    assert prepared.report["atoms_after"] > prepared.report["atoms_before"]
    assert any(" HD" in l or " HN" in l or l[12:14].strip().startswith("H")
               for l in prepared.pdb.read_text().splitlines() if l.startswith("ATOM"))


def test_additives_and_far_water_removed(prepared):
    removed = {c["residue"] for c in changes(prepared, "hetatm")}
    assert any("GOL" in r for r in removed) and any("SO4" in r for r in removed)
    assert "HOH" not in prepared.pdb.read_text()


def test_ligand_extracted_not_in_receptor(prepared):
    assert prepared.ligands and prepared.ligands[0]["resname"] == "LIG"
    assert "LIG" not in prepared.pdb.read_text()
    assert open(prepared.ligands[0]["path"]).read().count("HETATM") == 7


def test_zinc_in_site_kept_with_warning(prepared):
    assert any(c["residue"].endswith("ZN401") and "kept" in c["change"] for c in changes(prepared, "hetatm"))
    assert any("metal" in w for w in prepared.report["warnings"])
    assert "ZN" in prepared.pdbqt.read_text()


def test_protonation_recorded_with_tool(prepared):
    prot = changes(prepared, "protonation")
    assert prot and all(c["tool"] == "PDB2PQR+PROPKA" for c in prot)
    assert any("HIS" in c["change"] for c in prot)       # His -> HID/HIE/HIP
    assert prepared.report["net_charge_estimate"] == prepared.report["net_charge_estimate"]


def test_tools_versions_recorded(prepared):
    names = {t["tool"] for t in prepared.report["tools"]}
    assert {"pdbfixer", "pdb2pqr", "propka", "meeko"} <= names


def test_water_near_ligand_kept_on_request(decorated, tmp_path):
    r = prepare_receptor(decorated, tmp_path, ReceptorPrepOptions(reference_ligand="LIG", water_keep_radius=5.0))
    kept = [c for c in changes(r, "waters") if c["change"] == "kept"]
    assert len(kept) == 1 and "HOH201" in kept[0]["residue"]


def test_missing_atoms_rebuilt(tmp_path):
    st = peptide_structure()
    res = st[0][0][4]                      # truncate a side chain
    for n in [a.name for a in res if a.name not in ("N", "CA", "C", "O", "CB")]:
        res.remove_atom(n, "*")
    p = tmp_path / "t.pdb"; st.write_pdb(str(p))
    r = prepare_receptor(p, tmp_path / "o")
    assert any("added" in c["change"] for c in changes(r, "missing_atoms"))


def test_altloc_highest_occupancy(tmp_path):
    st = peptide_structure()
    res = st[0][0][3]
    a = res["CB"][0]; a.altloc, a.occ = "A", 0.3
    b = copy.copy(a); b.altloc, b.occ = "B", 0.7
    b.pos = b.pos.__class__(a.pos.x + 0.2, a.pos.y, a.pos.z); res.add_atom(b)
    p = tmp_path / "t.pdb"; st.write_pdb(str(p))
    r = prepare_receptor(p, tmp_path / "o")
    assert any("kept altloc B" in c["change"] for c in changes(r, "altlocs"))


def test_chain_break_warns(tmp_path):
    st = peptide_structure()
    del st[0][0][6]                         # remove a residue -> gap
    p = tmp_path / "t.pdb"; st.write_pdb(str(p))
    r = prepare_receptor(p, tmp_path / "o")
    assert any("chain break" in w for w in r.report["warnings"])


def test_flexible_residue_split(decorated, tmp_path):
    r = prepare_receptor(decorated, tmp_path, ReceptorPrepOptions(reference_ligand="LIG", flex_residues=["A:6"]))
    assert r.flex_pdbqt and r.flex_pdbqt.exists() and "BEGIN_RES" in r.flex_pdbqt.read_text()


def test_alphafold_plddt_warning_and_trim(tmp_path):
    st = peptide_structure()
    for i, res in enumerate(st[0][0]):
        for a in res:
            a.b_iso = 40.0 if i < 2 else 90.0
    p = tmp_path / "af.pdb"; st.write_pdb(str(p))
    r = prepare_receptor(p, tmp_path / "o", ReceptorPrepOptions(alphafold=True, trim_low_plddt_termini=True))
    assert len(changes(r, "alphafold")) == 2


def test_bad_inputs(tmp_path):
    (tmp_path / "bad.pdb").write_text("garbage\n")
    with pytest.raises(InvalidInput):
        prepare_receptor(tmp_path / "bad.pdb", tmp_path / "o")
    p = tmp_path / "t.pdb"; peptide_structure().write_pdb(str(p))
    with pytest.raises(InvalidInput):
        prepare_receptor(p, tmp_path / "o2", ReceptorPrepOptions(reference_ligand="XXX"))


def test_deterministic(decorated, tmp_path):
    o = ReceptorPrepOptions(reference_ligand="LIG")
    a = prepare_receptor(decorated, tmp_path / "a", o).pdbqt.read_text()
    b = prepare_receptor(decorated, tmp_path / "b", o).pdbqt.read_text()
    assert a == b


def test_minimize_restrains_heavy_atoms(tmp_path):
    p = tmp_path / "t.pdb"; peptide_structure().write_pdb(str(p))
    base = prepare_receptor(p, tmp_path / "a")
    mini = prepare_receptor(p, tmp_path / "b", ReceptorPrepOptions(minimize=True))
    assert not any("minimisation skipped" in w for w in mini.report["warnings"]), mini.report["warnings"]
    assert any(c["stage"] == "minimization" for c in mini.report["changes"])
    import gemmi
    a = gemmi.read_pdb_string(base.pdb.read_text())[0]
    b = gemmi.read_pdb_string(mini.pdb.read_text())[0]
    shifts = [x.pos.dist(y.pos) for ra, rb in zip(a[0], b[0]) for x, y in zip(ra, rb) if x.element.name != "H"]
    assert max(shifts) < 0.5          # heavy atoms stay near crystal coordinates


def test_mmcif_input(tmp_path):
    st = peptide_structure()
    st.make_mmcif_document().write_file(str(tmp_path / "t.cif"))
    r = prepare_receptor(tmp_path / "t.cif", tmp_path / "o")
    assert r.pdbqt.exists()
