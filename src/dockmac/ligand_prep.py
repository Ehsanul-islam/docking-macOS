"""Ligand preparation: bond orders, hydrogens, 3D coordinates, PDBQT.

Input can be a SMILES string, an SDF/MOL file, or a hetero group cut out of a
PDB file (which needs a SMILES template for its bond orders).
"""

from __future__ import annotations

import json
import urllib.error
import urllib.request
from pathlib import Path

import gemmi
from rdkit import Chem
from rdkit.Chem import AllChem


def from_smiles(smiles: str, optimize: bool = True) -> Chem.Mol:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"invalid SMILES: {smiles!r}")
    mol = Chem.AddHs(mol)
    if AllChem.EmbedMolecule(mol, randomSeed=42) != 0:
        raise ValueError(f"could not generate 3D coordinates for {smiles!r}")
    if optimize:
        AllChem.MMFFOptimizeMolecule(mol, maxIters=500)
    return mol


def from_file(path: str | Path) -> Chem.Mol:
    """Read the first molecule of an SDF/MOL/MOL2 file and add hydrogens."""
    path = Path(path)
    ext = path.suffix.lower()
    if ext == ".sdf":
        mol = next(iter(Chem.SDMolSupplier(str(path), removeHs=False)), None)
    elif ext == ".mol":
        mol = Chem.MolFromMolFile(str(path), removeHs=False)
    elif ext == ".mol2":
        mol = Chem.MolFromMol2File(str(path), removeHs=False)
    else:
        raise ValueError(f"unsupported ligand file type: {ext}")
    if mol is None:
        raise ValueError(f"could not read ligand from {path}")
    return _add_hs_3d(mol)


def fetch_template_smiles(resname: str, timeout: float = 15.0) -> str:
    """Look up the official SMILES of a PDB ligand code on RCSB (needs internet)."""
    url = f"https://data.rcsb.org/rest/v1/core/chemcomp/{resname.upper()}"
    try:
        with urllib.request.urlopen(url, timeout=timeout) as resp:
            data = json.load(resp)
    except (urllib.error.URLError, TimeoutError) as exc:
        raise ValueError(f"could not fetch template for {resname!r}: {exc}") from exc
    smiles = data.get("rcsb_chem_comp_descriptor", {}).get("smiles")
    if not smiles:
        raise ValueError(f"RCSB has no SMILES for {resname!r}")
    return smiles


def from_structure(ligand: gemmi.Structure, smiles_template: str) -> Chem.Mol:
    """Build a molecule from a ligand cut out of a PDB (see extract_ligand).

    A PDB ligand has no hydrogens or bond orders, so they cannot be reliably
    guessed from heavy-atom coordinates. The SMILES template (see
    ``fetch_template_smiles``) supplies bond orders and charges; hydrogens are
    then added on the crystal coordinates.
    """
    template = Chem.MolFromSmiles(smiles_template)
    if template is None:
        raise ValueError(f"invalid SMILES template: {smiles_template!r}")
    raw = Chem.MolFromPDBBlock(
        ligand.make_pdb_string(), removeHs=True, sanitize=False, proximityBonding=True
    )
    if raw is None:
        raise ValueError("could not parse ligand coordinates")
    try:
        mol = AllChem.AssignBondOrdersFromTemplate(template, raw)
    except ValueError as exc:
        raise ValueError("ligand does not match the SMILES template") from exc
    Chem.SanitizeMol(mol)
    if _flat_smiles(mol) != _flat_smiles(template):
        raise ValueError("ligand does not match the SMILES template")
    return _add_hs_3d(mol)


def _flat_smiles(mol: Chem.Mol) -> str:
    return Chem.MolToSmiles(Chem.RemoveHs(mol), isomericSmiles=False)


def _add_hs_3d(mol: Chem.Mol) -> Chem.Mol:
    """Add explicit hydrogens, placing them on the existing 3D coordinates."""
    has_3d = mol.GetNumConformers() > 0 and mol.GetConformer().Is3D()
    mol = Chem.AddHs(mol, addCoords=has_3d)
    if not has_3d:
        if AllChem.EmbedMolecule(mol, randomSeed=42) != 0:
            raise ValueError("could not generate 3D coordinates")
        AllChem.MMFFOptimizeMolecule(mol, maxIters=500)
    return mol


def to_pdbqt(mol: Chem.Mol) -> str:
    """Return AutoDock PDBQT text (Gasteiger charges, rotatable bonds)."""
    from meeko import MoleculePreparation, PDBQTWriterLegacy

    setups = MoleculePreparation().prepare(mol)
    pdbqt, ok, err = PDBQTWriterLegacy.write_string(setups[0])
    if not ok:
        raise ValueError(f"PDBQT conversion failed: {err}")
    return pdbqt


def write_sdf(mol: Chem.Mol, path: str | Path) -> None:
    with Chem.SDWriter(str(path)) as w:
        w.write(mol)
