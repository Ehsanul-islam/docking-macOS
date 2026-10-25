"""File formats and conversion (plan §7).

Everything converts through the canonical model (§4.5): gemmi Structure for macromolecules,
RDKit Mol for small molecules. Lossy conversions return warnings instead of failing silently (§7.2).
"""
from __future__ import annotations

from pathlib import Path

import gemmi
from rdkit import Chem
from rdkit.Chem import AllChem, rdDetermineBonds

from .errors import InvalidInput, NeedsTemplate, UnsupportedFeature
from .ligand_templates import is_untrusted, mark_untrusted

PDB_ATOM_LIMIT = 99_999
MIN_POLYMER_RESIDUES = 8        # fewer than this in a PDB file => treated as a small molecule
MACRO_EXT = {".pdb", ".ent", ".cif", ".mmcif"}
SMALL_EXT = {".sdf", ".mol", ".mol2", ".smi", ".smiles", ".xyz"}

# Sybyl atom typing for the MOL2 writer. RDKit has no MOL2 writer [VERIFY alternatives].
_ELEM_ONLY = {"H", "F", "Cl", "Br", "I", "Li", "Na", "K", "Ca", "Zn", "Mg", "Fe", "Cu", "Mn", "Co", "Ni"}


def sniff(path: str | Path) -> str:
    """Classify a file by content, not extension (drag-and-drop routing, plan §9.4).
    Returns: receptor | ligand | receptor_pdbqt | ligand_pdbqt | library | unknown."""
    path = Path(path)
    ext = path.suffix.lower()
    text = path.read_text(errors="replace")
    if ext == ".pdbqt":
        return "ligand_pdbqt" if ("ROOT" in text and "ENDROOT" in text) else "receptor_pdbqt"
    if ext in (".sdf", ".mol2"):
        n = text.count("$$$$") if ext == ".sdf" else text.count("@<TRIPOS>MOLECULE")
        return "library" if n > 1 else "ligand"
    if ext in (".smi", ".smiles"):
        return "library" if len([l for l in text.splitlines() if l.strip()]) > 1 else "ligand"
    if ext in (".mol", ".xyz"):
        return "ligand"
    if ext in MACRO_EXT:
        try:
            st = gemmi.read_structure(str(path))
        except Exception:  # noqa: BLE001
            return "unknown"
        n_res = 0
        for ch in st[0]:
            for r in ch:
                info = gemmi.find_tabulated_residue(r.name)
                if info is not None and (info.is_amino_acid() or info.is_nucleic_acid()):
                    n_res += 1
        return "receptor" if n_res >= MIN_POLYMER_RESIDUES else "ligand"
    return "unknown"


# ---------------------------------------------------------------- small molecules
def _sybyl_type(atom: Chem.Atom) -> str:
    el = atom.GetSymbol()
    if el in _ELEM_ONLY:
        return el
    arom = atom.GetIsAromatic()
    deg = atom.GetDegree() + atom.GetTotalNumHs()
    bonds = [b.GetBondType() for b in atom.GetBonds()]
    has_dbl = Chem.BondType.DOUBLE in bonds
    has_tpl = Chem.BondType.TRIPLE in bonds
    if el == "C":
        if arom:
            return "C.ar"
        return "C.1" if has_tpl or bonds.count(Chem.BondType.DOUBLE) == 2 else "C.2" if has_dbl else "C.3"
    if el == "N":
        if arom:
            return "N.ar"
        if atom.GetFormalCharge() == 1 and deg == 4 and not has_dbl:
            return "N.4"
        if has_tpl:
            return "N.1"
        if has_dbl:
            return "N.2"
        # amide / conjugated amine
        for nb in atom.GetNeighbors():
            if nb.GetSymbol() == "C" and any(b.GetBondType() == Chem.BondType.DOUBLE and
                                              b.GetOtherAtom(nb).GetSymbol() in ("O", "S")
                                              for b in nb.GetBonds()):
                return "N.am"
        return "N.pl3" if any(nb.GetIsAromatic() for nb in atom.GetNeighbors()) else "N.3"
    if el == "O":
        if has_dbl:
            return "O.2"
        if atom.GetFormalCharge() == -1 and any(
                nb.GetSymbol() in ("C", "S", "P") and sum(1 for x in nb.GetNeighbors() if x.GetSymbol() == "O") >= 2
                for nb in atom.GetNeighbors()):
            return "O.co2"
        return "O.3"
    if el == "S":
        n_o = sum(1 for nb in atom.GetNeighbors() if nb.GetSymbol() == "O")
        return "S.o2" if n_o >= 2 else "S.o" if n_o == 1 else "S.2" if has_dbl else "S.3"
    if el == "P":
        return "P.3"
    return el


def write_mol2(mol: Chem.Mol, path: str | Path, name: str | None = None) -> None:
    """Minimal Tripos MOL2 writer: first conformer, Gasteiger partial charges."""
    m = Chem.Mol(mol)
    if m.GetNumConformers() == 0:
        raise InvalidInput("MOL2 export needs 3D coordinates")
    Chem.Kekulize(m, clearAromaticFlags=False)
    AllChem.ComputeGasteigerCharges(m)
    conf = m.GetConformer()
    nm = name or (m.GetProp("_Name") if m.HasProp("_Name") and m.GetProp("_Name") else "MOL")
    lines = ["@<TRIPOS>MOLECULE", nm, f"{m.GetNumAtoms():5d}{m.GetNumBonds():5d}    0    0    0",
             "SMALL", "GASTEIGER", "", "@<TRIPOS>ATOM"]
    for a in m.GetAtoms():
        p = conf.GetAtomPosition(a.GetIdx())
        q = float(a.GetProp("_GasteigerCharge"))
        lines.append(f"{a.GetIdx() + 1:7d} {a.GetSymbol() + str(a.GetIdx() + 1):<8s}"
                     f"{p.x:10.4f}{p.y:10.4f}{p.z:10.4f} {_sybyl_type(a):<7s}1 {'LIG':<7s}{q:9.4f}")
    lines.append("@<TRIPOS>BOND")
    for i, b in enumerate(m.GetBonds(), 1):
        t = {Chem.BondType.SINGLE: "1", Chem.BondType.DOUBLE: "2", Chem.BondType.TRIPLE: "3",
             Chem.BondType.AROMATIC: "ar"}.get(b.GetBondType(), "1")
        if b.GetIsAromatic():
            t = "ar"
        lines.append(f"{i:6d}{b.GetBeginAtomIdx() + 1:6d}{b.GetEndAtomIdx() + 1:6d} {t}")
    Path(path).write_text("\n".join(lines) + "\n")


def read_molecules(path: str | Path, charge: int = 0) -> tuple[list[Chem.Mol], list[str]]:
    """Return (mols, warnings)."""
    path = Path(path)
    ext, warns = path.suffix.lower(), []
    if ext in (".sdf", ".mol"):
        mols = [m for m in Chem.SDMolSupplier(str(path), removeHs=False) if m is not None]
    elif ext == ".mol2":
        mols = []
        blocks = path.read_text().split("@<TRIPOS>MOLECULE")[1:]
        for b in blocks:
            m = Chem.MolFromMol2Block("@<TRIPOS>MOLECULE" + b, removeHs=False)
            if m is not None:
                mols.append(m)
    elif ext in (".smi", ".smiles"):
        mols = []
        for line in path.read_text().splitlines():
            if line.strip():
                m = Chem.MolFromSmiles(line.split()[0])
                if m is not None:
                    parts = line.split(None, 1)
                    if len(parts) > 1:
                        m.SetProp("_Name", parts[1].strip())
                    mols.append(m)
        warns.append("SMILES carry no 3D coordinates; 3D is generated during ligand preparation")
    elif ext == ".xyz":
        m = Chem.MolFromXYZFile(str(path))
        if m is None:
            raise InvalidInput("cannot read XYZ")
        try:
            rdDetermineBonds.DetermineBonds(m, charge=charge)
        except Exception as e:  # noqa: BLE001
            raise InvalidInput(f"bond perception failed for XYZ (charge {charge}): {e}") from e
        warns.append("XYZ has no bonds: bond orders and charges were perceived from geometry")
        mols = [m]
    elif ext == ".pdb":
        m = Chem.MolFromPDBFile(str(path), removeHs=False, proximityBonding=True)
        mols = [mark_untrusted(m)] if m else []
        warns.append("a PDB file has no bond orders or hydrogens: this molecule is a heavy-atom skeleton and "
                     "cannot be used for docking or written to SDF/MOL2/SMILES until bond orders are "
                     "established (use ligand.from_structure with a SMILES)")
    elif ext == ".pdbqt":
        mols, w = read_pdbqt_poses(path)
        warns += w
    else:
        raise UnsupportedFeature(f"unsupported small-molecule input {ext}")
    if not mols:
        raise InvalidInput(f"no molecules could be read from {path.name}")
    return mols, warns


def read_pdbqt_poses(path: str | Path) -> tuple[list[Chem.Mol], list[str]]:
    """PDBQT -> RDKit through Meeko's index map so bond orders are exact (plan §7.2 rule 3).
    Returns one Mol per pose (conformers split out), with `vina_score` set when present."""
    from meeko import PDBQTMolecule, RDKitMolCreate
    text = Path(path).read_text()
    if "REMARK SMILES" not in text:
        raise UnsupportedFeature(
            "PDBQT lacks Meeko's SMILES mapping, so bond orders cannot be recovered exactly; "
            "re-prepare the ligand with this app or supply the source SDF")
    pm = PDBQTMolecule(text, is_dlg=False, skip_typing=True)
    base = RDKitMolCreate.from_pdbqt_mol(pm)[0]
    scores = (pm._pose_data.get("free_energies") or [])
    mols = []
    for i, conf in enumerate(base.GetConformers()):
        m = Chem.Mol(base)
        m.RemoveAllConformers()
        c = Chem.Conformer(conf)
        c.SetId(0)
        m.AddConformer(c)
        m.SetIntProp("pose", i + 1)
        if i < len(scores):
            m.SetDoubleProp("vina_score", float(scores[i]))
        mols.append(m)
    return mols, []


def write_molecules(mols: list[Chem.Mol], path: str | Path) -> list[str]:
    path = Path(path)
    ext, warns = path.suffix.lower(), []
    if ext in (".sdf", ".mol", ".mol2", ".smi", ".smiles") and any(is_untrusted(m) for m in mols):
        raise NeedsTemplate("this molecule came from a PDB-style file and its bond orders are unknown; writing "
                            f"{ext} would record guessed chemistry. Establish bond orders first "
                            "(ligand.from_structure with a SMILES).", "", "smiles")
    if ext in (".sdf", ".mol"):
        if ext == ".mol" and len(mols) > 1:
            raise InvalidInput(".mol holds one molecule; use .sdf")
        w = Chem.SDWriter(str(path))
        for m in mols:
            w.write(m)
        w.close()
    elif ext == ".mol2":
        if len(mols) > 1:
            raise InvalidInput("multi-molecule MOL2 export is not supported; use SDF")
        write_mol2(mols[0], path)
        warns.append("MOL2 charges are Gasteiger partial charges; formal charges are encoded by Sybyl types only")
    elif ext in (".smi", ".smiles"):
        path.write_text("".join(
            f"{Chem.MolToSmiles(Chem.RemoveHs(m))} {m.GetProp('_Name') if m.HasProp('_Name') else ''}".rstrip() + "\n"
            for m in mols))
        warns.append("SMILES export drops 3D coordinates")
    elif ext == ".xyz":
        path.write_text(Chem.MolToXYZBlock(mols[0]))
        warns.append("XYZ export drops bonds, charges and stereochemistry")
    elif ext == ".pdb":
        Chem.MolToPDBFile(mols[0], str(path))
        warns.append("PDB export of a small molecule drops bond orders")
    else:
        raise UnsupportedFeature(f"unsupported small-molecule output {ext}")
    return warns


# ---------------------------------------------------------------- macromolecules
def read_structure(path: str | Path) -> tuple[gemmi.Structure, list[str]]:
    path = Path(path)
    warns = []
    if path.suffix.lower() == ".pdbqt":
        lines = []
        for l in path.read_text().splitlines():
            if l.startswith(("ATOM", "HETATM")):
                elem = l[77:79].strip() if len(l) >= 78 else ""
                elem = {"HD": "H", "OA": "O", "NA": "N", "SA": "S", "A": "C"}.get(elem, elem)
                lines.append(l[:54].ljust(76) + elem.rjust(2))
            elif l.startswith(("TER", "END")):
                lines.append(l)
        st = gemmi.read_pdb_string("\n".join(lines) + "\n")
        warns.append("PDBQT import: atom types collapsed to elements; partial charges, rotatable-bond "
                     "tree and bond orders are lost")
    else:
        try:
            st = gemmi.read_structure(str(path))
        except Exception as e:  # noqa: BLE001
            raise InvalidInput(f"cannot read {path.name}: {e}") from e
    return st, warns


def write_structure(st: gemmi.Structure, path: str | Path) -> list[str]:
    path = Path(path)
    ext, warns = path.suffix.lower(), []
    if ext in (".cif", ".mmcif"):
        st.make_mmcif_document().write_file(str(path))
    elif ext in (".pdb", ".ent"):
        n = sum(len(r) for ch in st[0] for r in ch)
        if n > PDB_ATOM_LIMIT:
            raise InvalidInput(f"{n} atoms exceed the PDB format limit ({PDB_ATOM_LIMIT}); export mmCIF")
        st.write_pdb(str(path))
    else:
        raise UnsupportedFeature(f"unsupported structure output {ext}")
    return warns


def convert(src: str | Path, dst: str | Path, charge: int = 0) -> dict:
    """Convert between supported formats. Returns {'warnings': [...], 'artifacts': [dst]}."""
    src, dst = Path(src), Path(dst)
    s_ext, d_ext = src.suffix.lower(), dst.suffix.lower()
    if s_ext == d_ext and s_ext not in (".pdb", ".sdf"):
        raise InvalidInput("source and destination formats are identical")
    kind = sniff(src)
    if kind in ("receptor", "receptor_pdbqt"):
        if d_ext == ".pdbqt":
            raise UnsupportedFeature("use receptor preparation to produce receptor PDBQT; direct conversion "
                                     "would skip hydrogens and charges")
        st, w1 = read_structure(src)
        w2 = write_structure(st, dst)
        return {"warnings": w1 + w2, "artifacts": [str(dst)]}
    if d_ext == ".pdbqt":
        raise UnsupportedFeature("use ligand preparation to produce ligand PDBQT")
    mols, w1 = read_molecules(src, charge)
    w2 = write_molecules(mols, dst)
    return {"warnings": w1 + w2, "artifacts": [str(dst)]}


def poses_to_sdf(pdbqt: str | Path, sdf: str | Path) -> dict:
    """Docked poses -> bond-order-correct multi-pose SDF with Vina scores."""
    mols, _ = read_pdbqt_poses(pdbqt)
    write_molecules(mols, sdf)
    return {"poses": len(mols), "artifacts": [str(sdf)]}
