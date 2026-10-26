"""Interaction profile geometry and the 2D renderer."""
import numpy as np
import pytest
from rdkit import Chem
from rdkit.Chem import AllChem

from dockchem import interactions, render2d


def ligand(smiles, seed=3):
    m = Chem.AddHs(Chem.MolFromSmiles(smiles))
    AllChem.EmbedMolecule(m, randomSeed=seed)
    return m


def pdb_atoms(rows):
    """rows: (serial, name, resn, resi, (x,y,z), element)"""
    return "\n".join(f"ATOM  {i:5d} {n:<4s} {r:>3s} A{k:4d}    {x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00          {e:>2s}"
                     for i, n, r, k, (x, y, z), e in rows) + "\nEND\n"


def pos(m, idx):
    return np.array(m.GetConformer().GetAtomPosition(idx))


def carbonyl_o(m):
    return next(a.GetIdx() for a in m.GetAtoms() if a.GetSymbol() == "O")


def prof(rec_rows, lig, **kw):
    return interactions.profile(pdb_atoms(rec_rows), lig, "T", "L", **kw)


def types(p):
    return sorted(i["type"] for i in p["interactions"])


def test_hydrogen_bond_with_a_linear_serine_donor():
    m = ligand("CC(=O)C")
    o = pos(m, carbonyl_o(m))
    u = np.array([1.0, 0, 0])
    p = prof([(1, "OG", "SER", 10, tuple(o + 2.8 * u), "O"), (2, "HG", "SER", 10, tuple(o + 1.85 * u), "H"),
              (3, "CB", "SER", 10, tuple(o + 3.9 * u + [0, 0.8, 0]), "C")], m)
    hb = [i for i in p["interactions"] if i["type"] == "Hydrogen bond"]
    assert len(hb) == 1 and hb[0]["residue"] == "Ser10" and hb[0]["ligand_role"] == "acceptor"
    assert hb[0]["dist"] == pytest.approx(2.8, abs=0.05)


def test_no_hydrogen_bond_beyond_3_5_angstrom_or_with_a_bent_donor():
    m = ligand("CC(=O)C")
    o = pos(m, carbonyl_o(m))
    far = prof([(1, "OG", "SER", 10, tuple(o + [3.8, 0, 0]), "O"), (2, "HG", "SER", 10, tuple(o + [2.85, 0, 0]), "H")], m)
    assert not [i for i in far["interactions"] if i["type"] == "Hydrogen bond"]
    # hydrogen pointing away from the acceptor (angle ~ 60 deg): not a hydrogen bond
    bent = prof([(1, "ND2", "ASN", 7, tuple(o + [3.0, 0, 0]), "N"), (2, "HD21", "ASN", 7, tuple(o + [3.5, 0.9, 0]), "H")], m)
    assert not [i for i in bent["interactions"] if i["type"] == "Hydrogen bond"]


def test_hydrophobic_contact_leucine_methyl():
    m = ligand("CCCC")
    c = next(a.GetIdx() for a in m.GetAtoms() if a.GetSymbol() == "C")
    p0 = pos(m, c)
    rows = [(1, "CB", "LEU", 5, tuple(p0 + [4.4, 3, 0]), "C"), (2, "CG", "LEU", 5, tuple(p0 + [4.4, 4.5, 0]), "C"),
            (3, "CD1", "LEU", 5, tuple(p0 + [3.5, 1.0, 0.5]), "C"), (4, "CD2", "LEU", 5, tuple(p0 + [5.5, 5.0, 0]), "C"),
            (5, "HA", "LEU", 5, tuple(p0 + [9, 9, 9]), "H")]
    p = prof(rows, m)
    assert "Hydrophobic" in types(p) and p["interactions"][0]["residue"] == "Leu5"


def benzene_ring(center, z, rot=0):
    ang = np.radians(np.arange(6) * 60 + rot)
    return [tuple(np.array(center) + [1.4 * np.cos(a), 1.4 * np.sin(a), z]) for a in ang]


def test_parallel_pi_stacking_and_no_stacking_when_far_or_offset():
    m = ligand("c1ccccc1")
    ring = [list(r) for r in m.GetRingInfo().AtomRings()][0]
    xyz = np.array([pos(m, i) for i in ring])
    c = xyz.mean(0)
    u, s, vt = np.linalg.svd(xyz - c)
    n = vt[2]
    a1 = np.cross(n, [1, 0, 0]); a1 = a1 / np.linalg.norm(a1)
    a2 = np.cross(n, a1)
    def phe(centre):
        return [(i + 1, nm, "PHE", 20, tuple(centre + 1.4 * (np.cos(np.radians(60 * i)) * a1 + np.sin(np.radians(60 * i)) * a2)), "C")
                for i, nm in enumerate(["CG", "CD1", "CE1", "CZ", "CE2", "CD2"])]
    far_h = (99, "HA", "PHE", 20, (60.0, 60.0, 60.0), "H")             # a hydrogen far away: the receptor is protonated
    def phe(centre, _p=phe):
        return _p(centre) + [far_h]
    stacked = prof(phe(c + 3.7 * n), m)
    assert "π-stacking (P)" in types(stacked)
    assert "π-stacking (P)" not in types(prof(phe(c + 7.0 * n), m))                    # too far
    assert "π-stacking (P)" not in types(prof(phe(c + 3.7 * n + 4.0 * a1), m))          # offset > 2 A


def test_salt_bridge_carboxylate_to_lysine():
    m = ligand("CC(=O)[O-]")
    ox = [a.GetIdx() for a in m.GetAtoms() if a.GetSymbol() == "O"]
    p0 = pos(m, ox[0])
    p = prof([(1, "NZ", "LYS", 30, tuple(p0 + [2.9, 0, 0]), "N"), (2, "HZ1", "LYS", 30, tuple(p0 + [1.9, 0, 0]), "H")], m)
    assert "Salt bridge" in types(p)


def test_van_der_waals_contact_only_when_nothing_specific():
    m = ligand("CCO")
    c = pos(m, 0)
    p = prof([(1, "CB", "SER", 40, tuple(c + [3.4, 0, 0]), "C"), (2, "OG", "SER", 40, tuple(c + [3.4, 1.4, 0]), "O"),
              (3, "HG", "SER", 40, tuple(c + [3.4, 2.4, 0.5]), "H")], m)
    assert [v["residue"] for v in p["vdw_contacts"]] == ["Ser40"]
    assert interactions.counts(p)["van der Waals"] == 1


def test_receptor_without_hydrogens_is_refused():
    with pytest.raises(ValueError, match="no hydrogens"):
        prof([(1, "CB", "SER", 40, (3.0, 0, 0), "C")], ligand("CCO"))


@pytest.fixture(scope="module")
def sample_profile():
    m = ligand("CC(=O)C")
    o = pos(m, carbonyl_o(m))
    return prof([(1, "OG", "SER", 10, tuple(o + [2.8, 0, 0]), "O"), (2, "HG", "SER", 10, tuple(o + [1.85, 0, 0]), "H"),
                 (3, "CB", "SER", 10, tuple(o + [3.9, 0.8, 0]), "C"), (4, "CB", "ALA", 11, tuple(o + [0.5, 3.6, 0]), "C")], m,
                affinity=-7.31)


@pytest.mark.parametrize("style", ["reference", "mono", "dark"])
def test_every_style_renders_a_png(tmp_path, sample_profile, style):
    out = tmp_path / f"{style}.png"
    n = render2d.draw(sample_profile, out, style=style, dpi=60)
    assert n >= 2 and out.stat().st_size > 5000 and out.read_bytes()[:4] == b"\x89PNG"


def test_render_handles_no_interactions_and_unknown_style(tmp_path):
    empty = dict(target_name="", compound="x", affinity=None, ligand_molblock=Chem.MolToMolBlock(ligand("CCO")),
                 interactions=[], vdw_contacts=[])
    assert render2d.draw(empty, tmp_path / "e.png", dpi=50) == 0
    with pytest.raises(KeyError):
        render2d.draw(empty, tmp_path / "x.png", style="nope")
