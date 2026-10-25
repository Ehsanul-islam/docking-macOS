import itertools

import gemmi
import pytest
from rdkit import Chem
from rdkit.Chem import AllChem

from dockchem import formats
from dockchem.errors import InvalidInput, UnsupportedFeature
from dockchem.ligprep import prepare_ligand

from fixtures import peptide_structure

# Neutral + charged + aromatic + stereo + amide + sulfonamide + carboxylate
SMILES = ["CC(=O)Nc1ccc(O)cc1", "C[C@H](N)C(=O)[O-]", "CS(=O)(=O)Nc1ccccc1", "C[N+](C)(C)CC#N", "O=C(O)c1ccncc1"]


def mol3d(smi):
    m = Chem.AddHs(Chem.MolFromSmiles(smi))
    assert AllChem.EmbedMolecule(m, randomSeed=5) == 0
    m.SetProp("_Name", "m")
    return m


def canon(m):
    return Chem.MolToSmiles(Chem.RemoveHs(m))


@pytest.mark.parametrize("smi", SMILES)
def test_sdf_roundtrip_exact(tmp_path, smi):
    m = mol3d(smi)
    formats.write_molecules([m], tmp_path / "a.sdf")
    back, _ = formats.read_molecules(tmp_path / "a.sdf")
    assert canon(back[0]) == canon(m)
    assert back[0].GetNumAtoms() == m.GetNumAtoms()
    assert abs(back[0].GetConformer().GetPositions() - m.GetConformer().GetPositions()).max() < 1e-3


@pytest.mark.parametrize("smi", SMILES)
def test_mol2_roundtrip_topology(tmp_path, smi):
    m = mol3d(smi)
    formats.write_molecules([m], tmp_path / "a.mol2")
    back, _ = formats.read_molecules(tmp_path / "a.mol2")
    assert back[0].GetNumAtoms() == m.GetNumAtoms() and back[0].GetNumBonds() == m.GetNumBonds()
    assert sum(a.GetFormalCharge() for a in back[0].GetAtoms()) == sum(a.GetFormalCharge() for a in m.GetAtoms())
    assert abs(back[0].GetConformer().GetPositions() - m.GetConformer().GetPositions()).max() < 1e-3


def test_mol2_needs_3d(tmp_path):
    with pytest.raises(InvalidInput):
        formats.write_molecules([Chem.MolFromSmiles("CCO")], tmp_path / "a.mol2")


def test_smiles_roundtrip_warns(tmp_path):
    m = mol3d("C[C@H](N)C(=O)[O-]")
    w = formats.write_molecules([m], tmp_path / "a.smi")
    back, w2 = formats.read_molecules(tmp_path / "a.smi")
    assert canon(back[0]) == canon(m) and w and w2


def test_xyz_roundtrip_warns_and_perceives_bonds(tmp_path):
    m = mol3d("CC(=O)Nc1ccc(O)cc1")
    w = formats.write_molecules([m], tmp_path / "a.xyz")
    back, w2 = formats.read_molecules(tmp_path / "a.xyz")
    assert any("drops bonds" in x for x in w) and any("perceived" in x for x in w2)
    assert back[0].GetNumAtoms() == m.GetNumAtoms()
    assert canon(back[0]) == canon(m)


def test_all_small_pairs_via_convert(tmp_path):
    m = mol3d("CC(=O)Nc1ccc(O)cc1")
    formats.write_molecules([m], tmp_path / "src.sdf")
    for ext in (".mol2", ".xyz", ".smi", ".mol"):
        r = formats.convert(tmp_path / "src.sdf", tmp_path / f"o{ext}")
        assert (tmp_path / f"o{ext}").exists()
        assert isinstance(r["warnings"], list)


def test_pdbqt_to_sdf_keeps_bond_orders_and_scores(tmp_path):
    p = prepare_ligand("CC(=O)Nc1ccc(O)cc1", "x")[0]
    pq = tmp_path / "out.pdbqt"
    pq.write_text("MODEL 1\nREMARK VINA RESULT:    -6.500      0.000      0.000\n" + p.pdbqt + "ENDMDL\n"
                  "MODEL 2\nREMARK VINA RESULT:    -6.100      1.200      2.000\n" + p.pdbqt + "ENDMDL\n")
    r = formats.poses_to_sdf(pq, tmp_path / "poses.sdf")
    assert r["poses"] == 2
    mols = list(Chem.SDMolSupplier(str(tmp_path / "poses.sdf"), removeHs=False))
    assert canon(mols[0]) == "CC(=O)Nc1ccc(O)cc1"
    assert [m.GetDoubleProp("vina_score") for m in mols] == [-6.5, -6.1]


def test_pdbqt_without_mapping_refused(tmp_path):
    (tmp_path / "x.pdbqt").write_text("ROOT\nATOM      1  C   UNL     1       0.0 0.0 0.0  1.00  0.00    +0.000 C\nENDROOT\nTORSDOF 0\n")
    with pytest.raises(UnsupportedFeature):
        formats.read_pdbqt_poses(tmp_path / "x.pdbqt")


def test_structure_pdb_cif_roundtrip(tmp_path):
    st = peptide_structure()
    st.write_pdb(str(tmp_path / "a.pdb"))
    formats.convert(tmp_path / "a.pdb", tmp_path / "a.cif")
    formats.convert(tmp_path / "a.cif", tmp_path / "b.pdb")
    a, b = gemmi.read_structure(str(tmp_path / "a.pdb")), gemmi.read_structure(str(tmp_path / "b.pdb"))
    pa = [(x.name, x.pos.x, x.pos.y, x.pos.z) for r in a[0][0] for x in r]
    pb = [(x.name, x.pos.x, x.pos.y, x.pos.z) for r in b[0][0] for x in r]
    assert len(pa) == len(pb) and all(abs(u[1] - v[1]) < 1e-3 for u, v in zip(pa, pb))


def test_pdb_atom_limit(tmp_path):
    st = peptide_structure()
    big = gemmi.Structure()
    m = gemmi.Model("1")
    for i in range(70):                      # 70 x ~ 15 residues x ~8 atoms > 99,999 is too heavy; fake via loop
        pass
    n_needed = formats.PDB_ATOM_LIMIT // sum(len(r) for r in st[0][0]) + 2
    for i in range(n_needed):
        ch = gemmi.Chain(f"C{i}")
        for r in st[0][0]:
            ch.add_residue(r)
        m.add_chain(ch)
    big.add_model(m)
    with pytest.raises(InvalidInput):
        formats.write_structure(big, tmp_path / "big.pdb")
    formats.write_structure(big, tmp_path / "big.cif")        # mmCIF has no such limit


def test_pdbqt_receptor_import_is_lossy(tmp_path):
    pq = tmp_path / "r.pdbqt"
    pq.write_text("ATOM      1  N   GLY A   1       1.000   2.000   3.000  1.00  0.00    -0.350 N \n"
                  "ATOM      2  HN  GLY A   1       1.500   2.500   3.500  1.00  0.00    +0.160 HD\n")
    st, w = formats.read_structure(pq)
    assert st[0][0][0][1].element.name == "H" and any("lost" in x for x in w)


def test_receptor_pdbqt_export_refused(tmp_path):
    peptide_structure().write_pdb(str(tmp_path / "a.pdb"))
    with pytest.raises(UnsupportedFeature):
        formats.convert(tmp_path / "a.pdb", tmp_path / "a.pdbqt")


def test_sniff_routes_by_content(tmp_path):
    peptide_structure().write_pdb(str(tmp_path / "prot.pdb"))
    assert formats.sniff(tmp_path / "prot.pdb") == "receptor"
    formats.write_molecules([mol3d("CCO"), mol3d("CCN")], tmp_path / "lib.sdf")
    assert formats.sniff(tmp_path / "lib.sdf") == "library"
    Chem.MolToPDBFile(mol3d("CC(=O)Nc1ccc(O)cc1"), str(tmp_path / "lig.pdb"))
    assert formats.sniff(tmp_path / "lig.pdb") == "ligand"
    p = prepare_ligand("CCO", "e")[0]
    (tmp_path / "l.pdbqt").write_text(p.pdbqt)
    assert formats.sniff(tmp_path / "l.pdbqt") == "ligand_pdbqt"
