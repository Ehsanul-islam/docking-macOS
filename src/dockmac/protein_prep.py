"""Protein preparation, the DockPrep equivalent.

Steps: keep best alternate location, fill missing heavy atoms, convert
nonstandard residues, add hydrogens at a chosen pH. Waters, ions and ligands
are handled beforehand with ``structure.clean``.
"""

from __future__ import annotations

import io
from dataclasses import dataclass
from pathlib import Path

import gemmi
from openmm.app import PDBFile
from pdbfixer import PDBFixer


@dataclass
class PrepReport:
    missing_atoms: int
    missing_residues_added: int
    nonstandard_replaced: int
    hydrogens_added: int


def prepare_protein(
    st: gemmi.Structure,
    ph: float = 7.0,
    add_hydrogens: bool = True,
    fill_missing_residues: bool = False,
) -> tuple[gemmi.Structure, PrepReport]:
    """Return (prepared structure, report). Input is not modified.

    Missing residues are only rebuilt when ``fill_missing_residues`` is set,
    and then only internal gaps, never missing chain ends.
    """
    buf = io.StringIO(st.make_pdb_string())
    fixer = PDBFixer(pdbfile=buf)

    fixer.findMissingResidues()
    if fill_missing_residues:
        chains = list(fixer.topology.chains())
        for key in list(fixer.missingResidues):
            chain = chains[key[0]]
            if key[1] == 0 or key[1] >= len(list(chain.residues())):
                del fixer.missingResidues[key]  # terminal gap: leave alone
    else:
        fixer.missingResidues = {}
    n_res = sum(len(v) for v in fixer.missingResidues.values())

    fixer.findNonstandardResidues()
    n_nonstd = len(fixer.nonstandardResidues)
    fixer.replaceNonstandardResidues()

    fixer.removeHeterogens(keepWater=False)  # ligands are prepared separately
    fixer.findMissingAtoms()
    n_atoms = sum(len(v) for v in fixer.missingAtoms.values())
    fixer.addMissingAtoms()

    before = fixer.topology.getNumAtoms()
    if add_hydrogens:
        fixer.addMissingHydrogens(ph)
    n_h = fixer.topology.getNumAtoms() - before if add_hydrogens else 0

    out = io.StringIO()
    PDBFile.writeFile(fixer.topology, fixer.positions, out, keepIds=True)
    prepared = gemmi.read_pdb_string(out.getvalue())
    prepared.setup_entities()
    return prepared, PrepReport(n_atoms, n_res, n_nonstd, n_h)


def save(st: gemmi.Structure, path: str | Path) -> None:
    st.write_pdb(str(path))
