"""Bond orders for ligands taken from structure files (plan §6 stage 1, §8.1).

PDB/mmCIF coordinates carry no bond orders and usually no hydrogens, so a ligand cut out of a crystal
structure is only a heavy-atom skeleton. Guessing orders from that skeleton gives plausible-looking but
wrong molecules (imatinib came out saturated, biotin's acid became a gem-diol), so this module refuses to
guess. A ligand gets bond orders from exactly one of:

  1. a SMILES supplied by the user,
  2. the bundled dictionary (`data/ligand_dictionary.json`), whose entries were each validated against a
     real crystal ligand, or
  3. perception from geometry + net charge, only when the file already contains hydrogens.

Otherwise NeedsTemplate is raised, saying what to supply.
"""
from __future__ import annotations

import json
from collections import Counter
from dataclasses import dataclass, field
from pathlib import Path

import gemmi
from rdkit import Chem
from rdkit.Chem import AllChem, rdDetermineBonds

from .errors import InvalidInput, NeedsTemplate

DICT_PATH = Path(__file__).parent / "data" / "ligand_dictionary.json"
TRUST = "_bond_orders_trusted"          # mol property: "1" established, "0" unknown
SOURCE = "_bond_orders_source"

# Expected lengths (A) for multiple bonds, used only as a sanity audit of an assigned template.
DOUBLE_RANGES = {("C", "O"): (1.10, 1.34), ("C", "N"): (1.18, 1.42), ("C", "C"): (1.26, 1.46),
                 ("O", "P"): (1.38, 1.62), ("O", "S"): (1.36, 1.60), ("N", "O"): (1.14, 1.34)}
AUDIT_SLACK = 0.12


@dataclass
class ResolvedLigand:
    mol: Chem.Mol                       # explicit hydrogens, 3D, bond orders established
    source: str                         # user_smiles | dictionary:<RES> | perceived_with_charge
    warnings: list[str] = field(default_factory=list)
    audit: list[str] = field(default_factory=list)


def load_dictionary() -> dict:
    return json.loads(DICT_PATH.read_text())["entries"]


def mark_untrusted(mol: Chem.Mol) -> Chem.Mol:
    mol.SetProp(TRUST, "0")
    return mol


def is_untrusted(mol: Chem.Mol) -> bool:
    return mol.HasProp(TRUST) and mol.GetProp(TRUST) == "0"


def residue_to_mol(path: str | Path, resname: str, chain: str | None = None,
                   seq: int | None = None) -> tuple[Chem.Mol, gemmi.Residue, list[str]]:
    """Cut one residue out of a structure file. Heavy-atom skeleton, single bonds, flagged untrusted."""
    try:
        st = gemmi.read_structure(str(path))
    except Exception as e:  # noqa: BLE001
        raise InvalidInput(f"cannot read {Path(path).name}: {e}") from e
    found = [(ch.name, r) for ch in st[0] for r in ch
             if r.name == resname and (chain is None or ch.name == chain)
             and (seq is None or r.seqid.num == seq)]
    if not found:
        raise InvalidInput(f"{resname} not found in {Path(path).name}")
    warns = []
    if len(found) > 1:
        c, r = found[0]
        warns.append(f"{len(found)} copies of {resname}; using {c}:{r.seqid.num}")
    cname, res = found[0]
    # keep the first altloc only
    s2 = gemmi.Structure(); m = gemmi.Model("1"); c = gemmi.Chain(cname)
    r2 = gemmi.Residue(); r2.name, r2.seqid, r2.het_flag = res.name, res.seqid, "H"
    seen = set()
    for a in res:
        if a.name in seen:
            continue
        seen.add(a.name)
        r2.add_atom(a)
    c.add_residue(r2); m.add_chain(c); s2.add_model(m)
    mol = Chem.MolFromPDBBlock(s2.make_pdb_string(), proximityBonding=True, removeHs=False, sanitize=False)
    if mol is None or mol.GetNumAtoms() == 0:
        raise InvalidInput(f"cannot build a molecule from {resname}")
    return mark_untrusted(mol), res, warns


def _heavy_formula(mol: Chem.Mol) -> Counter:
    return Counter(a.GetSymbol() for a in mol.GetAtoms() if a.GetAtomicNum() > 1)


def audit_bond_lengths(mol: Chem.Mol) -> list[str]:
    """Flag double bonds whose length is implausible for their atom pair (sign of a misassigned order)."""
    conf = mol.GetConformer()
    out = []
    for b in mol.GetBonds():
        if b.GetBondType() != Chem.BondType.DOUBLE:
            continue
        a, c = b.GetBeginAtom(), b.GetEndAtom()
        key = tuple(sorted((a.GetSymbol(), c.GetSymbol())))
        rng = DOUBLE_RANGES.get(key)
        if not rng:
            continue
        d = (conf.GetAtomPosition(a.GetIdx()) - conf.GetAtomPosition(c.GetIdx())).Length()
        if not rng[0] - AUDIT_SLACK <= d <= rng[1] + AUDIT_SLACK:
            out.append(f"{a.GetSymbol()}{a.GetIdx() + 1}={c.GetSymbol()}{c.GetIdx() + 1} is {d:.2f} A")
    return out


def _finish(mol: Chem.Mol, source: str) -> Chem.Mol:
    """Stereo from the crystal geometry (carbon centres only: 3D 'stereo' on phosphates is an artefact)."""
    Chem.AssignStereochemistryFrom3D(mol)
    for a in mol.GetAtoms():
        if a.GetAtomicNum() != 6 and a.GetChiralTag() != Chem.ChiralType.CHI_UNSPECIFIED:
            a.SetChiralTag(Chem.ChiralType.CHI_UNSPECIFIED)
    mol.SetProp(TRUST, "1")
    mol.SetProp(SOURCE, source)
    return mol


def apply_template(mol: Chem.Mol, smiles: str, resname: str, source: str) -> Chem.Mol:
    template = Chem.MolFromSmiles(smiles)
    if template is None:
        raise InvalidInput(f"cannot parse SMILES: {smiles}")
    have, want = _heavy_formula(mol), _heavy_formula(template)
    if have != want:
        diff = {k: (have.get(k, 0), want.get(k, 0)) for k in set(have) | set(want) if have.get(k, 0) != want.get(k, 0)}
        raise NeedsTemplate(
            f"{resname}: the structure has {sum(have.values())} heavy atoms but the template has "
            f"{sum(want.values())} (element counts structure/template: {diff}). The ligand may be incomplete in "
            "the file, covalently modified, or the template may be for a different molecule.",
            resname, "smiles")
    mol.UpdatePropertyCache(strict=False)
    try:
        assigned = AllChem.AssignBondOrdersFromTemplate(template, mol)
    except Exception as e:  # noqa: BLE001
        raise NeedsTemplate(
            f"{resname}: the template does not match the structure's connectivity ({e}). Check the SMILES.",
            resname, "smiles") from e
    # Hydrogens with coordinates, so explicit-H geometry matches the crystal heavy atoms.
    Chem.SanitizeMol(assigned)
    out = Chem.AddHs(assigned, addCoords=True)
    return _finish(out, source)


def _perceive_with_charge(mol: Chem.Mol, charge: int, resname: str) -> Chem.Mol:
    if not any(a.GetAtomicNum() == 1 for a in mol.GetAtoms()):
        raise NeedsTemplate(
            f"{resname}: the file has no hydrogens, so bond orders cannot be derived from a net charge "
            "alone. Supply a SMILES for this ligand.", resname, "smiles")
    m = Chem.Mol(mol)
    try:
        rdDetermineBonds.DetermineBonds(m, charge=charge)
    except Exception as e:  # noqa: BLE001
        raise NeedsTemplate(
            f"{resname}: bond perception failed for net charge {charge} ({e}). Check the charge or supply a SMILES.",
            resname, "smiles") from e
    return _finish(m, "perceived_with_charge")


def ligand_from_structure(path: str | Path, resname: str, chain: str | None = None, seq: int | None = None,
                          smiles: str | None = None, charge: int | None = None) -> ResolvedLigand:
    """Establish bond orders for a ligand in a structure file, or raise NeedsTemplate."""
    mol, res, warns = residue_to_mol(path, resname, chain, seq)
    if smiles:
        out, source = apply_template(mol, smiles, resname, "user_smiles"), "user_smiles"
    else:
        entry = load_dictionary().get(resname)
        if entry:
            out = apply_template(mol, entry["smiles"], resname, f"dictionary:{resname}")
            source = f"dictionary:{resname}"
            warns.append(f"bond orders and protonation state for {resname} come from the bundled dictionary "
                         f"({entry.get('protonation', 'neutral form')})")
        elif charge is not None:
            out, source = _perceive_with_charge(mol, charge, resname), "perceived_with_charge"
        else:
            raise NeedsTemplate(
                f"{resname} is not in the bundled ligand dictionary, so its bond orders are unknown. "
                "Supply a SMILES (or, if the file contains hydrogens, the net charge).", resname, "smiles")
    audit = audit_bond_lengths(out)
    if audit:
        warns.append(f"{len(audit)} double bond(s) have implausible lengths; check the template: "
                     + "; ".join(audit[:4]))
    out.SetProp("_Name", f"{resname}")
    return ResolvedLigand(out, source, warns, audit)
