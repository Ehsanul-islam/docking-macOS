"""Synthetic receptor fixtures. Real PDB entries can't be fetched in CI sandboxes, so we build a
peptide with RDKit and decorate it the way crystal structures are decorated."""
import gemmi
from rdkit import Chem
from rdkit.Chem import AllChem

SEQ = "GSHDEKRCNQFYTAL"          # covers His/Asp/Glu/Lys/Arg/Cys/Asn/Gln


def peptide_structure(seq=SEQ, seed=2) -> gemmi.Structure:
    m = Chem.AddHs(Chem.MolFromSequence(seq))
    ps = AllChem.ETKDGv3()
    ps.randomSeed, ps.useRandomCoords, ps.maxIterations = seed, True, 2000
    assert AllChem.EmbedMolecule(m, ps) == 0
    AllChem.MMFFOptimizeMolecule(m, maxIters=5000)
    st = gemmi.read_pdb_string(Chem.MolToPDBBlock(Chem.RemoveHs(m), flavor=0))
    st.setup_entities()
    return st


def _add_het(chain, resname, atoms, seqid, het="H"):
    res = gemmi.Residue()
    res.name, res.seqid, res.het_flag = resname, gemmi.SeqId(seqid, " "), het
    for name, el, pos in atoms:
        a = gemmi.Atom()
        a.name, a.element, a.pos, a.occ, a.b_iso = name, gemmi.Element(el), gemmi.Position(*pos), 1.0, 20.0
        res.add_atom(a)
    chain.add_residue(res)


def decorated_structure() -> gemmi.Structure:
    """Peptide + ligand 'LIG' near residue centre + far water/glycerol/sulfate/zinc."""
    st = peptide_structure()
    ch = st[0][0]
    ca = [r["CA"][0].pos for r in ch if r.find_atom("CA", "*")]
    c = ca[len(ca) // 2]
    cen = [sum(p.x for p in ca) / len(ca), sum(p.y for p in ca) / len(ca), sum(p.z for p in ca) / len(ca)]
    d = [c.x - cen[0], c.y - cen[1], c.z - cen[2]]
    n = sum(v * v for v in d) ** 0.5
    x, y, z = (c.x + 4.5 * d[0] / n, c.y + 4.5 * d[1] / n, c.z + 4.5 * d[2] / n)   # just outside the surface
    _add_het(ch, "LIG", [("C1", "C", (x, y, z)), ("C2", "C", (x + 1.5, y, z)),
                         ("C3", "C", (x + 1.5, y + 1.5, z)), ("N1", "N", (x, y + 1.5, z)),
                         ("O1", "O", (x + 3.0, y, z)), ("C4", "C", (x, y, z + 1.5)),
                         ("C5", "C", (x + 3.0, y + 1.5, z))], 101)
    _add_het(ch, "HOH", [("O", "O", (x + 2, y + 2, z + 2))], 201)           # near ligand
    _add_het(ch, "HOH", [("O", "O", (c.x + 40, c.y, c.z))], 202)            # bulk water
    _add_het(ch, "GOL", [("C1", "C", (c.x + 30, c.y, c.z)), ("O1", "O", (c.x + 31.4, c.y, c.z))], 301)
    _add_het(ch, "SO4", [("S", "S", (c.x - 30, c.y, c.z)), ("O1", "O", (c.x - 31.4, c.y, c.z))], 302)
    _add_het(ch, "ZN", [("ZN", "Zn", (x - 2, y - 2, z - 2))], 401)          # in site
    return st
