"""Protein-ligand interaction profile for one docked pose (feeds the 2D interaction image, render2d.py).

Rules follow the published figure pipeline this app reproduces (every label comes from this profile, nothing is drawn
by hand):
  hydrogen bonds   heavy-atom D...A <= 3.5 A and D-H...A >= 120 deg; rotatable OH / NH3+ donors (hydrogens not
                   optimised by the engine): X-D...A >= 90 deg
  salt bridges     ligand carboxylate (assumed ionised at pH 7.4) to Lys NZ / Arg NE,NH1,NH2 <= 4.0 A
  water bridges    ligand polar atom - kept water <= 3.5 A - protein N/O <= 3.5 A
  hydrophobic      ligand C (only C/H neighbours) to protein C of a hydrophobic residue, 0.5-4.0 A
  pi-stacking      aromatic ring centres <= 5.5 A, parallel (<30 deg) or T-shaped (60-90 deg), offset <= 2.0 A
  pi-cation        aromatic ligand ring - Lys NZ / Arg guanidinium <= 6.0 A, offset <= 2.0 A
  van der Waals    residues with a heavy atom <= 4.0 A of the ligand and no specific interaction
The receptor must carry hydrogens (the app's prepared receptor does); the ligand is a 3D RDKit molecule.
"""
from __future__ import annotations

import numpy as np
from rdkit import Chem

THREE = {"ALA": "Ala", "ARG": "Arg", "ASN": "Asn", "ASP": "Asp", "CYS": "Cys", "GLN": "Gln", "GLU": "Glu", "GLY": "Gly",
         "HIS": "His", "HID": "His", "HIE": "His", "HIP": "His", "ILE": "Ile", "LEU": "Leu", "LYS": "Lys", "MET": "Met",
         "PHE": "Phe", "PRO": "Pro", "SER": "Ser", "THR": "Thr", "TRP": "Trp", "TYR": "Tyr", "VAL": "Val", "HOH": "HOH",
         "CYX": "Cys", "HEM": "Haem ", "ASH": "Asp", "GLH": "Glu", "LYN": "Lys", "CYM": "Cys"}
HYDROPHOBIC_RES = {"ALA", "VAL", "LEU", "ILE", "MET", "PHE", "TRP", "PRO", "TYR"}
RINGS = {"PHE": [("CG", "CD1", "CD2", "CE1", "CE2", "CZ")], "TYR": [("CG", "CD1", "CD2", "CE1", "CE2", "CZ")],
         "TRP": [("CG", "CD1", "CD2", "NE1", "CE2"), ("CD2", "CE2", "CE3", "CZ2", "CZ3", "CH2")],
         "HIS": [("CG", "ND1", "CD2", "CE1", "NE2")], "HID": [("CG", "ND1", "CD2", "CE1", "NE2")],
         "HIE": [("CG", "ND1", "CD2", "CE1", "NE2")], "HIP": [("CG", "ND1", "CD2", "CE1", "NE2")]}
ROT = {("SER", "OG"), ("THR", "OG1"), ("TYR", "OH"), ("LYS", "NZ"), ("LYN", "NZ")}


def _angle(a, b, c) -> float:
    v1, v2 = a - b, c - b
    return float(np.degrees(np.arccos(np.clip(v1 @ v2 / (np.linalg.norm(v1) * np.linalg.norm(v2) + 1e-9), -1, 1))))


def _normal(pts):
    c = pts.mean(0)
    u, s, vt = np.linalg.svd(pts - c)
    return c, vt[2]


def read_receptor(pdb_text: str, cofactors=("HEM",)) -> list[dict]:
    atoms = []
    for l in pdb_text.splitlines():
        if l.startswith("ATOM") or (l.startswith("HETATM") and l[17:20].strip() in cofactors + ("HOH",)):
            name = l[12:16].strip()
            atoms.append(dict(name=name, resn=l[17:20].strip(), resi=int(l[22:26]), chain=l[21],
                              el=(l[76:78].strip() or name[0]).capitalize() if l[76:78].strip() != "D" else "H",
                              xyz=np.array([float(l[30:38]), float(l[38:46]), float(l[46:54])])))
    return atoms


def profile(receptor_pdb: str, ligand: Chem.Mol, target_name: str = "", compound: str = "ligand",
            affinity: float | None = None, score_label: str = "AutoDock Vina binding affinity") -> dict:
    """Interaction profile dict (the schema render2d.draw() consumes)."""
    prot = read_receptor(receptor_pdb)
    if not any(a["el"] == "H" for a in prot):
        raise ValueError("the receptor has no hydrogens: hydrogen bonds cannot be assigned reliably")
    mol = Chem.Mol(ligand)
    if not any(a.GetAtomicNum() == 1 for a in mol.GetAtoms()):
        mol = Chem.AddHs(mol, addCoords=True)
    lxyz = np.array([list(mol.GetConformer().GetAtomPosition(i)) for i in range(mol.GetNumAtoms())])
    chains = {a["chain"] for a in prot if a["resn"] != "HOH"}
    multichain = len(chains) > 1

    def resid(rt, rn, ch):
        return f"{THREE.get(rt, rt.title())}{rn}" + (f" ({ch})" if multichain else "")

    ints: list[dict] = []

    def add(kind, rt, rn, ch, lig_xyz, prot_xyz, dist, extra=None):
        ints.append(dict(type=kind, residue=resid(rt, rn, ch), restype=rt, resnr=int(rn), chain=ch,
                         lig=[float(x) for x in lig_xyz], prot=[float(x) for x in prot_xyz],
                         dist=round(float(dist), 2), **(extra or {})))

    P = np.array([a["xyz"] for a in prot])
    heavy = [i for i, a in enumerate(prot) if a["el"] != "H"]
    # --- hydrogens attached to their parent heavy atom
    Hs: dict[int, list] = {}
    for i, a in enumerate(prot):
        if a["el"] == "H":
            d = np.linalg.norm(P - a["xyz"], axis=1)
            d[i] = 9
            j = int(np.argmin(d))
            if d[j] < 1.25:
                Hs.setdefault(j, []).append(a["xyz"])
    by_res: dict = {}
    for i in heavy:
        a = prot[i]
        by_res.setdefault((a["chain"], a["resi"]), {})[a["name"]] = i
    heavy_xyz = P[heavy]

    def bonded(i, j) -> bool:
        return np.linalg.norm(P[i] - P[j]) < 1.9

    # --- hydrophobic contacts
    lig_c = [a.GetIdx() for a in mol.GetAtoms() if a.GetSymbol() == "C"
             and all(n.GetSymbol() in ("C", "H") for n in a.GetNeighbors())]
    prot_c = []
    for i in heavy:
        a = prot[i]
        if a["el"] == "C" and a["resn"] in HYDROPHOBIC_RES and a["name"] not in ("C", "CA"):
            nbrs = [j for j in heavy if j != i and prot[j]["chain"] == a["chain"] and prot[j]["resi"] == a["resi"] and bonded(i, j)]
            if all(prot[j]["el"] == "C" for j in nbrs):
                prot_c.append(i)
    if lig_c and prot_c:
        D = np.linalg.norm(lxyz[lig_c][:, None, :] - P[prot_c][None, :, :], axis=2)
        for ki, pi_ in zip(*np.where((D > 0.5) & (D <= 4.0))):
            a = prot[prot_c[pi_]]
            add("Hydrophobic", a["resn"], a["resi"], a["chain"], lxyz[lig_c[ki]], P[prot_c[pi_]], D[ki, pi_])

    # --- ring terms: pi-stacking and pi-cation, genuinely aromatic ligand rings only
    lig_rings = [list(r) for r in mol.GetRingInfo().AtomRings() if all(mol.GetAtomWithIdx(i).GetIsAromatic() for i in r)]
    for r in lig_rings:
        lc, ln = _normal(lxyz[r])
        for (ch, resi), names in by_res.items():
            rn = prot[next(iter(names.values()))]["resn"]
            for ring in RINGS.get(rn, []):
                if not all(n in names for n in ring):
                    continue
                pc, pn = _normal(P[[names[n] for n in ring]])
                dist = float(np.linalg.norm(lc - pc))
                if dist > 5.5:
                    continue
                ang = float(np.degrees(np.arccos(np.clip(abs(ln @ pn), 0, 1))))
                off1 = float(np.linalg.norm((pc - lc) - ((pc - lc) @ ln) * ln))
                off2 = float(np.linalg.norm((lc - pc) - ((lc - pc) @ pn) * pn))
                kind = None
                if ang < 30 and min(off1, off2) <= 2.0:
                    kind = "π-stacking (P)"
                elif 60 <= ang <= 90 and min(off1, off2) <= 2.0:
                    kind = "π-stacking (T)"
                if kind:
                    add(kind, rn, resi, ch, lc, pc, dist, dict(lig_ring=[list(map(float, lxyz[i])) for i in r]))
        for i in heavy:
            a = prot[i]
            cat = None
            if a["resn"] in ("LYS", "LYN") and a["name"] == "NZ" and a["resn"] == "LYS":
                cat = P[i]
            elif a["resn"] == "ARG" and a["name"] == "CZ":
                cat = P[i]
            if cat is None:
                continue
            dist = float(np.linalg.norm(lc - cat))
            off = float(np.linalg.norm((cat - lc) - ((cat - lc) @ ln) * ln))
            if dist <= 6.0 and off <= 2.0:
                add("π-cation", a["resn"], a["resi"], a["chain"], lc, cat, dist, dict(lig_ring=[list(map(float, lxyz[k])) for k in r]))

    # --- hydrogen bonds, salt bridges, water bridges
    def p_donor(i): return prot[i]["el"] in ("N", "O") and i in Hs
    def p_acceptor(i):
        a = prot[i]
        if a["el"] == "O":
            return True
        return a["el"] == "N" and a["resn"] in ("HIS", "HID", "HIE") and a["name"] in ("ND1", "NE2") and i not in Hs
    lig_H = {a.GetIdx(): [lxyz[n.GetIdx()] for n in a.GetNeighbors() if n.GetAtomicNum() == 1] for a in mol.GetAtoms()}
    def l_donor(k): return mol.GetAtomWithIdx(k).GetSymbol() in ("N", "O") and len(lig_H[k]) > 0
    def l_acceptor(k):
        a = mol.GetAtomWithIdx(k)
        if a.GetSymbol() == "O":
            return True
        return a.GetSymbol() == "N" and not lig_H[k] and a.GetDegree() < 4

    def hb(dxyz, hlist, axyz, ante=None):
        dist = float(np.linalg.norm(dxyz - axyz))
        if dist > 3.5 or not hlist:
            return None
        if ante is not None:
            ang = _angle(ante, dxyz, axyz)
            return (dist, ang) if ang >= 90 else None
        best = max(_angle(dxyz, h, axyz) for h in hlist)
        return (dist, best) if best >= 120 else None

    def p_ante(i):
        a = prot[i]
        if (a["resn"], a["name"]) not in ROT:
            return None
        d = np.linalg.norm(P - P[i], axis=1)
        cand = [j for j in heavy if j != i and prot[j]["resi"] == a["resi"] and prot[j]["chain"] == a["chain"] and d[j] < 1.7]
        return P[cand[0]] if cand else None

    def l_ante(k):
        a = mol.GetAtomWithIdx(k)
        hn = [n.GetIdx() for n in a.GetNeighbors() if n.GetAtomicNum() > 1]
        return lxyz[hn[0]] if a.GetSymbol() == "O" and len(hn) == 1 else None

    lpolar = [a.GetIdx() for a in mol.GetAtoms() if a.GetSymbol() in ("N", "O")]
    near_prot = [i for i in heavy if prot[i]["el"] in ("N", "O")
                 and lpolar and np.linalg.norm(lxyz[lpolar] - P[i], axis=1).min() <= 3.5]
    for k in lpolar:
        for i in near_prot:
            a = prot[i]
            if a["resn"] == "HOH":
                continue
            if l_donor(k) and p_acceptor(i):
                r = hb(lxyz[k], lig_H[k], P[i], l_ante(k))
                if r:
                    add("Hydrogen bond", a["resn"], a["resi"], a["chain"], lxyz[k], P[i], r[0],
                        dict(ligand_role="donor", prot_atom=a["name"], angle=round(r[1], 1)))
            if l_acceptor(k) and p_donor(i):
                r = hb(P[i], Hs[i], lxyz[k], p_ante(i))
                if r:
                    add("Hydrogen bond", a["resn"], a["resi"], a["chain"], lxyz[k], P[i], r[0],
                        dict(ligand_role="acceptor", prot_atom=a["name"], angle=round(r[1], 1)))
    patt = Chem.MolFromSmarts("[CX3](=O)[OX2H1,OX1-]")
    for m in mol.GetSubstructMatches(patt):
        ox = [m[1], m[2]]
        for i in heavy:
            a = prot[i]
            if (a["resn"] in ("LYS", "LYN") and a["name"] == "NZ") or (a["resn"] == "ARG" and a["name"] in ("NE", "NH1", "NH2")):
                d = min(np.linalg.norm(lxyz[o] - P[i]) for o in ox)
                if d <= 4.0:
                    add("Salt bridge", a["resn"], a["resi"], a["chain"], lxyz[ox].mean(0), P[i], d,
                        dict(prot_atom=a["name"], note="carboxylate assumed ionised at pH 7.4"))
    for w in [i for i in heavy if prot[i]["resn"] == "HOH"]:
        dl = [t for t in ((float(np.linalg.norm(lxyz[k] - P[w])), k) for k in lpolar) if t[0] <= 3.5]
        if not dl:
            continue
        dlw, k = min(dl)
        for i in [j for j in heavy if prot[j]["el"] in ("N", "O") and prot[j]["resn"] != "HOH" and np.linalg.norm(P[j] - P[w]) <= 3.5]:
            a = prot[i]
            add("Water bridge", a["resn"], a["resi"], a["chain"], lxyz[k], P[i], dlw,
                dict(water=[float(x) for x in P[w]], water_name=f"HOH{prot[w]['resi']}", prot_atom=a["name"],
                     water_prot_dist=round(float(np.linalg.norm(P[i] - P[w])), 2)))
    seen, uniq = set(), []
    for i in ints:
        key = (i["type"], i["residue"], tuple(np.round(i["lig"], 2)), tuple(np.round(i["prot"], 2)))
        if key not in seen:
            seen.add(key)
            uniq.append(i)
    ints = uniq

    # --- van der Waals contacts (no specific interaction)
    lig_heavy = np.array([lxyz[a.GetIdx()] for a in mol.GetAtoms() if a.GetAtomicNum() > 1])
    near: dict = {}
    for i in heavy:
        a = prot[i]
        if a["resn"] == "HOH":
            continue
        dv = np.linalg.norm(lig_heavy - P[i], axis=1)
        d = float(dv.min())
        if d <= 4.0:
            k = (a["resn"], a["resi"], a["chain"])
            if d < near.get(k, (9, None))[0]:
                near[k] = (d, lig_heavy[int(np.argmin(dv))])
    specific = {(i["restype"], i["resnr"], i["chain"]) for i in ints}
    contacts = [dict(residue=resid(rt, rn, ch), restype=rt, resnr=rn, chain=ch, dist=round(float(d), 2),
                     lig=[float(x) for x in la])
                for (rt, rn, ch), (d, la) in sorted(near.items(), key=lambda t: t[0][1]) if (rt, rn, ch) not in specific]
    return dict(target_name=target_name, compound=compound, affinity=affinity, score_label=score_label,
                ligand_molblock=Chem.MolToMolBlock(mol), interactions=ints, vdw_contacts=contacts)


def counts(prof: dict) -> dict:
    out: dict = {}
    for i in prof["interactions"]:
        out[i["type"]] = out.get(i["type"], 0) + 1
    out["van der Waals"] = len(prof["vdw_contacts"])
    return out
