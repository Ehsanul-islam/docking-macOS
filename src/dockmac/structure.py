"""Load a structure, describe its contents, and delete chains / hetero groups.

Everything works on a copy; the input structure is never modified.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import gemmi

WATER_NAMES = {"HOH", "WAT", "DOD", "H2O", "TIP", "TIP3", "SOL"}


@dataclass
class HeteroGroup:
    chain: str
    resname: str
    seqid: str
    n_atoms: int
    kind: str  # "water", "ion" or "ligand"


@dataclass
class ChainInfo:
    name: str
    n_polymer_residues: int
    hetero: list[HeteroGroup] = field(default_factory=list)


def load(path: str | Path) -> gemmi.Structure:
    """Read a PDB/mmCIF file, keep the first model only."""
    st = gemmi.read_structure(str(path))
    st.setup_entities()
    while len(st) > 1:
        del st[len(st) - 1]
    return st


def classify(residue: gemmi.Residue) -> str:
    """Return "polymer", "water", "ion" or "ligand" for a residue."""
    if residue.name in WATER_NAMES:
        return "water"
    if residue.entity_type == gemmi.EntityType.Polymer:
        return "polymer"
    info = gemmi.find_tabulated_residue(residue.name)
    if info is not None and (info.is_amino_acid() or info.is_nucleic_acid()):
        return "polymer"
    if len(residue) == 1:
        return "ion"
    return "ligand"


def summarize(st: gemmi.Structure) -> list[ChainInfo]:
    """List every chain with its polymer size and hetero groups."""
    out = []
    for chain in st[0]:
        info = ChainInfo(chain.name, 0)
        for res in chain:
            kind = classify(res)
            if kind == "polymer":
                info.n_polymer_residues += 1
            else:
                info.hetero.append(
                    HeteroGroup(chain.name, res.name, str(res.seqid), len(res), kind)
                )
        out.append(info)
    return out


def _keep_best_altloc(st: gemmi.Structure) -> None:
    """Per residue keep the altloc with the highest total occupancy."""
    for chain in st[0]:
        for res in chain:
            totals: dict[str, float] = {}
            for atom in res:
                if atom.has_altloc():
                    totals[atom.altloc] = totals.get(atom.altloc, 0.0) + atom.occ
            if len(totals) < 2:
                continue
            best = max(sorted(totals), key=lambda a: totals[a])
            drop = [(a.name, a.altloc) for a in res if a.has_altloc() and a.altloc != best]
            for name, altloc in drop:
                res.remove_atom(name, altloc)
    st.remove_alternative_conformations()  # clears altloc ids, drops leftovers


def clean(
    st: gemmi.Structure,
    keep_chains: list[str] | None = None,
    remove_water: bool = True,
    remove_ions: bool = True,
    remove_ligands: bool = False,
    keep_ligands: list[str] | None = None,
    remove_hydrogens: bool = False,
    best_altloc: bool = True,
) -> gemmi.Structure:
    """Return a cleaned copy of ``st``.

    keep_chains: chain names to keep (None keeps all).
    keep_ligands: residue names that survive ``remove_ligands`` (e.g. ["ERL"]).
    """
    new = st.clone()
    keep_lig = {n.upper() for n in (keep_ligands or [])}
    model = new[0]

    if keep_chains is not None:
        wanted = set(keep_chains)
        for name in [c.name for c in model if c.name not in wanted]:
            model.remove_chain(name)

    for chain in model:
        drop = []
        for i, res in enumerate(chain):
            kind = classify(res)
            if (
                (kind == "water" and remove_water)
                or (kind == "ion" and remove_ions and res.name not in keep_lig)
                or (kind == "ligand" and remove_ligands and res.name not in keep_lig)
            ):
                drop.append(i)
        for i in reversed(drop):
            del chain[i]

    for name in [c.name for c in model if len(c) == 0]:
        model.remove_chain(name)

    if best_altloc:
        _keep_best_altloc(new)
    if remove_hydrogens:
        new.remove_hydrogens()
    new.setup_entities()
    return new


def extract_ligand(
    st: gemmi.Structure, resname: str, chain: str | None = None
) -> gemmi.Structure:
    """Return a new structure holding only the first ``resname`` hetero group."""
    for ch in st[0]:
        if chain is not None and ch.name != chain:
            continue
        for res in ch:
            if res.name == resname and classify(res) in ("ligand", "ion"):
                new = gemmi.Structure()
                new.cell, new.spacegroup_hm = st.cell, st.spacegroup_hm
                model = gemmi.Model("1")
                new_chain = gemmi.Chain(ch.name)
                new_chain.add_residue(res)
                model.add_chain(new_chain)
                new.add_model(model)
                return new
    raise ValueError(f"ligand {resname!r} not found")


def write_pdb(st: gemmi.Structure, path: str | Path) -> None:
    st.write_pdb(str(path))
