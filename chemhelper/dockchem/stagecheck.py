"""Output checks run after each preparation stage, so a failure is reported where it starts.

Each check raises ToolFailure naming the stage and the specific problems; none of them repairs anything.
"""
from __future__ import annotations

import math
import re
from pathlib import Path

import gemmi

from .errors import ToolFailure

AD_TYPES = {"C", "A", "N", "NA", "OA", "O", "S", "SA", "HD", "H", "P", "F", "Cl", "Br", "I", "Zn", "Mg", "Mn",
            "Ca", "Fe", "B", "Cu", "Ni", "Co", "Cd", "Hg"}
HEAVY_MOVE_LIMIT = 1.0      # A: PDB2PQR legitimately nudges a few side chains (<= ~0.9 A seen); more is a defect
PDBFIXER_MOVE_LIMIT = 0.05  # A: PDBFixer must not move atoms that were already there


def _atoms(text: str) -> dict:
    out = {}
    for l in text.splitlines():
        if l.startswith(("ATOM", "HETATM")):
            out[(l[21], l[22:27].strip(), l[12:16].strip())] = (
                float(l[30:38]), float(l[38:46]), float(l[46:54]), l[17:20].strip(), l[76:78].strip())
    return out


def _fail(stage: str, problems: list[str]):
    err = ToolFailure(f"{stage} produced an invalid result: " + "; ".join(problems[:5])
                      + (f" (+{len(problems) - 5} more)" if len(problems) > 5 else ""),
                      log="\n".join(problems))
    err.stage = stage
    raise err


def residue_set(text: str) -> dict:
    return {(c, n): r for (c, n, _), (_, _, _, r, _) in _atoms(text).items()}


def check_moved_atoms(stage: str, before: str, after: str, allowed_moved: set | None = None,
                      limit: float = HEAVY_MOVE_LIMIT) -> float:
    """Residues must survive and their input heavy atoms must stay put (except `allowed_moved` residue keys,
    e.g. deliberately flipped ones). Returns the largest displacement seen."""
    a, b = _atoms(before), _atoms(after)
    problems, worst = [], 0.0
    lost = set(residue_set(before)) - set(residue_set(after))
    if lost:
        problems.append("residues lost: " + ", ".join(f"{c}:{n}" for c, n in sorted(lost)))
    for k, (x, y, z, res, el) in a.items():
        if el == "H" or k[2].startswith("H") or (k[0], k[1]) in (allowed_moved or set()):
            continue
        q = b.get(k)
        if q is None:
            continue                      # removed atoms are reported by the stage that removes them
        d = math.dist((x, y, z), q[:3])
        worst = max(worst, d)
        if d > limit:
            problems.append(f"{k[0]}:{res}{k[1]}.{k[2]} moved {d:.2f} A")
    if problems:
        _fail(stage, problems)
    return worst


def check_hydrogens(stage: str, text: str):
    """Every amino-acid residue of a protonated model must carry hydrogens, and none may be non-finite."""
    problems = []
    per = {}
    for (c, n, name), (x, y, z, res, el) in _atoms(text).items():
        if not all(map(math.isfinite, (x, y, z))):
            problems.append(f"{c}:{res}{n}.{name} has non-finite coordinates")
        per.setdefault((c, n, res), [0, 0])[0 if el == "H" else 1] += 1
    for (c, n, res), (h, heavy) in per.items():
        info = gemmi.find_tabulated_residue(res)
        if info is not None and info.is_amino_acid() and h == 0:
            problems.append(f"{c}:{res}{n} has no hydrogens")
    if problems:
        _fail(stage, problems)


def check_pdbqt(prepared_pdb: str, pdbqt: Path, dropped: list[dict], flex: Path | None = None):
    """The PDBQT must hold every non-dropped residue's heavy atoms, with known finite atom types."""
    prepared = _atoms(prepared_pdb)
    gone = set()
    for d in dropped:
        chain, _, rest = d["residue"].partition(":")
        m = re.search(r"(-?\d+[A-Za-z]?)$", rest)
        if m:
            gone.add((chain, m.group(1)))
    out = _atoms_pdbqt(Path(pdbqt))
    if flex and Path(flex).exists():
        out.update(_atoms_pdbqt(Path(flex)))        # flexible side chains live in their own file
    problems = []
    for k, (x, y, z, res, ad) in out.items():
        if ad not in AD_TYPES:
            problems.append(f"{k[0]}:{res}{k[1]}.{k[2]} has unknown atom type {ad!r}")
        if not all(map(math.isfinite, (x, y, z))):
            problems.append(f"{k[0]}:{res}{k[1]}.{k[2]} has non-finite coordinates")
    heavy_in, heavy_out = {}, {}
    for k, (_, _, _, res, el) in prepared.items():
        if el != "H":
            heavy_in.setdefault((k[0], k[1]), (res, set()))[1].add(k[2])
    for k, (_, _, _, res, ad) in out.items():
        if ad not in ("HD", "H"):
            heavy_out.setdefault((k[0], k[1]), set()).add(k[2])
    for key, (res, names) in sorted(heavy_in.items()):
        if key in gone:
            continue
        got = heavy_out.get(key)
        if got is None:
            problems.append(f"{key[0]}:{res}{key[1]} is absent from the PDBQT without being reported as dropped")
        elif len(got) < len(names - {"OXT"}):          # names can differ (Meeko renames water atoms): compare counts
            problems.append(f"{key[0]}:{res}{key[1]} lost heavy atoms in the PDBQT ({len(got)} of {len(names)})")
    if not out:
        problems.append("the PDBQT is empty")
    if problems:
        _fail("pdbqt", problems)
    return len(out)


def _atoms_pdbqt(path: Path) -> dict:
    out = {}
    for l in path.read_text().splitlines():
        if l.startswith(("ATOM", "HETATM")):
            out[(l[21], l[22:27].strip(), l[12:16].strip())] = (
                float(l[30:38]), float(l[38:46]), float(l[46:54]), l[17:20].strip(), l[77:79].strip())
    return out
