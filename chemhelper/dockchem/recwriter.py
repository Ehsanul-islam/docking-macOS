"""Meeko-free rigid-receptor PDBQT writer (fallback for gap 5).

Meeko's receptor step perceives bonds from atom distances, so unusual geometry can make it fail. For a rigid
Vina receptor all that is needed is coordinates and one AutoDock atom type per atom, and for the standard
residues that type follows from (residue name, atom name). This writer looks the type up in a table
harvested from Meeko's own output (tools/receptor_types/harvest.py) and checked against Meeko in the tests,
so on a structure both can handle they agree atom for atom.

What it does not do, deliberately:
- guess a type: an atom with no table entry makes its residue "unsupported", which is reported and
  recorded as a dropped residue (status degraded), never typed by a heuristic;
- write partial charges: Vina ignores them, they are written as 0.000 (AutoDock4 scoring would need Meeko);
- write flexible side chains (needs a torsion tree): the caller must use Meeko for those.
"""
from __future__ import annotations

import json
from functools import lru_cache
from pathlib import Path

import gemmi

from .errors import InvalidInput

TABLE = Path(__file__).parent / "data" / "receptor_atom_types.json"
BOND_ALLOWANCE = 1.2    # Meeko bonds two atoms closer than this x the sum of covalent radii
SS_MAX = BOND_ALLOWANCE * 2 * 1.05      # A: SG-SG distance Meeko reads as a disulfide (S radius 1.05)
_METALS = {"ZN": "Zn", "MN": "Mn", "MG": "Mg", "CA": "Ca", "FE": "Fe"}


@lru_cache(maxsize=1)
def _table() -> dict:
    return json.loads(TABLE.read_text())


def supported_residues() -> set[str]:
    return set(_table()["types"])


def _parent_element(hpos, res) -> str:
    """Element of the heavy atom a hydrogen is attached to (nearest within 1.3 A), else ''."""
    best, el = 1.3, ""
    for a in res:
        if a.element.name != "H":
            d = hpos.dist(a.pos)
            if d < best:
                best, el = d, a.element.name
    return el


def write_pdbqt(pdb_text: str, out: Path) -> dict:
    """Write `out` and return {"atoms": n, "dropped": [{"residue", "reason"}], "merged_hydrogens": n}."""
    t = _table()
    types = t["types"]
    st = gemmi.read_pdb_string(pdb_text)
    if not len(st):
        raise InvalidInput("no model in the receptor")
    model = st[0]
    sg = [(ch.name, r.seqid, a.pos) for ch in model for r in ch for a in r if r.name in ("CYS", "CYX") and a.name == "SG"]

    def in_disulfide(chain, seqid, pos):
        return any((c, s) != (chain, seqid) and pos.dist(p) < SS_MAX for c, s, p in sg)

    lines, dropped, n_merged, serial = [], [], 0, 0
    for ch in model:
        for res in ch:
            label = f"{ch.name}:{res.name}{res.seqid.num}{res.seqid.icode.strip()}"
            table = types.get(res.name)
            if table is None and res.name.upper() in _METALS and len(res) == 1:
                table = {res.name.upper(): _METALS[res.name.upper()]}
            if table is None:
                dropped.append({"residue": label, "reason": f"no atom-type rules for {res.name}"})
                continue
            rows, missing = [], None
            for a in res:
                ad = table.get(a.name)
                if ad is None and a.name == "OXT" and a.element.name == "O":
                    ad = "OA"                                   # C-terminal carboxylate oxygen
                if ad is None and a.element.name == "H":
                    parent = _parent_element(a.pos, res)        # AutoDock rule: polar H kept (HD), others merged
                    if parent == "C":
                        n_merged += 1
                        continue
                    if parent in ("N", "O", "S"):
                        ad = "HD"
                if ad is None:
                    missing = a.name
                    break
                if res.name == "CYX" and a.name == "SG":        # two heavy neighbours (CB, SG) => acceptor
                    ad = "SA" if in_disulfide(ch.name, res.seqid, a.pos) else "S"
                rows.append((a, ad))
            if missing:
                dropped.append({"residue": label, "reason": f"no atom-type rule for {res.name}.{missing}"})
                continue
            for a, ad in rows:
                serial += 1
                name = a.name if len(a.name) == 4 else f" {a.name:<3s}"
                lines.append(f"ATOM  {serial % 100000:5d} {name} {res.name:>3s} {ch.name[:1]}{res.seqid.num:4d}"
                             f"{res.seqid.icode if res.seqid.icode.strip() else ' '}   "
                             f"{a.pos.x:8.3f}{a.pos.y:8.3f}{a.pos.z:8.3f}{a.occ:6.2f}{0.0:6.2f}    "
                             f"{0.0:+6.3f} {ad:<2s}")
    Path(out).write_text("\n".join(lines) + "\n")
    return {"atoms": serial, "dropped": dropped, "merged_hydrogens": n_merged}


# ---- cofactors, sugars and other hetero groups kept as part of the receptor ---------------------------------

def type_hetero_group(pdb_block: str, resname: str, chain: str, seq: int) -> tuple[list[tuple], list[str]]:
    """AutoDock-type one kept hetero group (FAD, NAG, NAP, ...): bond orders from the bundled wwPDB dictionary
    (ligand_templates), hydrogens added by RDKit, atom types by Meeko's own typing. Returns
    ([(name, (x, y, z), ad_type)], warnings). Raises NeedsTemplate / InvalidInput when bond orders or element
    support are missing; the caller then reports the group as left out. Metal complexes are refused."""
    import tempfile

    from rdkit import Chem
    from rdkit.Chem import rdForceFieldHelpers  # noqa: F401  (ensures RDKit force-field data is loaded)

    from .ligand_templates import ligand_from_structure
    from .ligprep import COMMON_ELEMENTS

    with tempfile.TemporaryDirectory() as td:
        p = Path(td) / "group.pdb"
        p.write_text(pdb_block)
        resolved = ligand_from_structure(p, resname, chain, seq)
    mol = resolved.mol
    odd = sorted({a.GetSymbol() for a in mol.GetAtoms()} - COMMON_ELEMENTS)
    if odd:
        raise InvalidInput(f"{resname} contains {', '.join(odd)}: no reliable AutoDock atom types")
    names = {a.GetIdx(): (a.GetPDBResidueInfo().GetName().strip() if a.GetPDBResidueInfo() else f"{a.GetSymbol()}{a.GetIdx()}")
             for a in mol.GetAtoms()}
    molh = Chem.AddHs(mol, addCoords=True)
    from meeko import MoleculePreparation
    try:
        setup = MoleculePreparation().prepare(molh)[0]
    except Exception:  # noqa: BLE001  Gasteiger fails on hypervalent centres; Vina ignores charges
        setup = MoleculePreparation(charge_model="zero").prepare(molh)[0]
    conf = molh.GetConformer()
    out, nh = [], 0
    for a in molh.GetAtoms():
        i = a.GetIdx()
        if setup.get_is_ignore(i):                  # non-polar hydrogens are merged into their carbon
            continue
        ad = setup.get_atom_type(i)
        if not ad:
            raise InvalidInput(f"{resname}: atom {names.get(i, i)} got no AutoDock type")
        if a.GetAtomicNum() == 1:
            nh += 1
            name = f"H{nh}"
        else:
            name = names.get(i, f"{a.GetSymbol()}{i}")
        p3 = conf.GetAtomPosition(i)
        out.append((name, (p3.x, p3.y, p3.z), ad))
    return out, list(resolved.warnings) + [f"hydrogens of {resname} were added from its dictionary protonation "
                                             f"state ({resolved.source}), not from the structure"]


def append_hetero(out: Path, chain: str, resname: str, seq: int, atoms: list[tuple], serial0: int) -> int:
    """Append typed atoms of one hetero group to a PDBQT; returns the next serial number."""
    lines, serial = [], serial0
    for name, (x, y, z), ad in atoms:
        serial += 1
        nm = name if len(name) == 4 else f" {name:<3s}"
        lines.append(f"ATOM  {serial % 100000:5d} {nm} {resname:>3s} {chain[:1]}{seq:4d}    "
                     f"{x:8.3f}{y:8.3f}{z:8.3f}{1.0:6.2f}{0.0:6.2f}    {0.0:+6.3f} {ad:<2s}")
    with open(out, "a") as fh:
        fh.write("\n".join(lines) + "\n")
    return serial
