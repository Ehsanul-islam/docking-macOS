"""Write a receptor PDBQT from a prepared (hydrogenated) protein.

AutoDock atom types follow the AutoDockTools convention: aromatic carbons are
``A``, acceptor oxygens ``OA``, sulfur ``SA``, unprotonated histidine
nitrogens ``NA``, polar hydrogens ``HD``. Non-polar hydrogens are merged into
their carbon. Partial charges come from the Amber ff14SB force field; if a
residue has no force-field template the charges fall back to zero (Vina does
not use them).
"""

from __future__ import annotations

import io
from dataclasses import dataclass
from pathlib import Path

import gemmi
from openmm import NonbondedForce
from openmm.app import ForceField, PDBFile
from openmm.unit import elementary_charge

_AROMATIC = {
    "PHE": {"CG", "CD1", "CD2", "CE1", "CE2", "CZ"},
    "TYR": {"CG", "CD1", "CD2", "CE1", "CE2", "CZ"},
    "TRP": {"CG", "CD1", "CD2", "CE2", "CE3", "CZ2", "CZ3", "CH2"},
    "HIS": {"CG", "CD2", "CE1"},
    "HID": {"CG", "CD2", "CE1"},
    "HIE": {"CG", "CD2", "CE1"},
    "HIP": {"CG", "CD2", "CE1"},
}
_HIS = {"HIS", "HID", "HIE", "HIP"}
_POLAR_PARENTS = {"N", "O", "S"}
_OTHER = {"ZN": "Zn", "MG": "Mg", "MN": "Mn", "FE": "Fe", "CA": "Ca", "CL": "Cl", "F": "F", "BR": "Br", "I": "I", "P": "P"}
_BOND = 1.3  # max H-parent distance in angstrom


@dataclass
class ReceptorReport:
    n_atoms: int
    charges_from_forcefield: bool


def _charges(st: gemmi.Structure) -> list[float] | None:
    """Per-atom charges in file order, or None if the force field can't type it."""
    pdb = PDBFile(io.StringIO(st.make_pdb_string()))
    try:
        system = ForceField("amber14-all.xml").createSystem(pdb.topology)
    except ValueError:
        return None
    nb = next(f for f in system.getForces() if isinstance(f, NonbondedForce))
    return [nb.getParticleParameters(i)[0].value_in_unit(elementary_charge) for i in range(system.getNumParticles())]


def _parent(res: gemmi.Residue, h: gemmi.Atom) -> gemmi.Atom | None:
    best, best_d = None, _BOND
    for a in res:
        if a.is_hydrogen():
            continue
        d = a.pos.dist(h.pos)
        if d < best_d:
            best, best_d = a, d
    return best


def to_pdbqt(st: gemmi.Structure) -> tuple[str, ReceptorReport]:
    """Return (pdbqt text, report) for a prepared protein structure."""
    charges = _charges(st)
    lines: list[str] = []
    serial = 0
    idx = 0  # index into the charge list (file order)
    for chain in st[0]:
        for res in chain:
            atoms = list(res)
            charge_of = {}
            for a in atoms:
                charge_of[id(a)] = charges[idx] if charges else 0.0
                idx += 1
            # parent bookkeeping: H bound to polar atoms stay (HD), others merge
            keep_h, merged = {}, {}
            has_h = set()
            for a in atoms:
                if a.is_hydrogen():
                    p = _parent(res, a)
                    if p is None:
                        continue
                    has_h.add(p.name)
                    if p.element.name in _POLAR_PARENTS:
                        keep_h[a.name] = p
                    else:
                        merged[p.name] = merged.get(p.name, 0.0) + charge_of[id(a)]
            for a in atoms:
                el = a.element.name
                if a.is_hydrogen():
                    if a.name not in keep_h:
                        continue
                    ad, q = "HD", charge_of[id(a)]
                else:
                    q = charge_of[id(a)] + merged.get(a.name, 0.0)
                    if el == "C":
                        ad = "A" if a.name in _AROMATIC.get(res.name, ()) else "C"
                    elif el == "N":
                        ad = "NA" if res.name in _HIS and a.name in ("ND1", "NE2") and a.name not in has_h else "N"
                    elif el == "O":
                        ad = "OA"
                    elif el == "S":
                        ad = "SA"
                    elif el.upper() in _OTHER:
                        ad = _OTHER[el.upper()]
                    else:
                        raise ValueError(f"no AutoDock type for element {el} in {res.name}")
                serial += 1
                name = a.name if len(a.name) >= 4 else f" {a.name}"
                lines.append(
                    f"ATOM  {serial:5d} {name:<4s} {res.name:>3s} {chain.name[:1]:1s}{res.seqid.num:4d}    "
                    f"{a.pos.x:8.3f}{a.pos.y:8.3f}{a.pos.z:8.3f}{1.0:6.2f}{0.0:6.2f}    {q:+6.3f} {ad:<2s}"
                )
    return "\n".join(lines) + "\nEND\n", ReceptorReport(serial, charges is not None)


def write(st: gemmi.Structure, path: str | Path) -> ReceptorReport:
    text, report = to_pdbqt(st)
    Path(path).write_text(text)
    return report
