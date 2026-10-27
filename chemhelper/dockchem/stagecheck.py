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


ANCHORS = {"N", "CA", "C", "O", "CB"}       # atoms a preparation stage must not relocate
ANCHOR_LIMIT = 0.5                          # A
BOND_LIMIT = 1.9                            # A: heavy-atom pairs closer than this in the input are bonded
BOND_CHANGE = 0.35                          # A: bonds must keep their length (rotations and flips do)


def check_moved_atoms(stage: str, before: str, after: str, allowed_moved: set | None = None,
                      limit: float = ANCHOR_LIMIT) -> float:
    """What a stage may not do: lose residues, relocate the backbone / CB, or stretch or squash bonds.
    Side chains legitimately rotate (flips, PDB2PQR's debumping) and equivalent atoms (OD1/OD2, ...) are
    relabelled, which is why positions of other atoms are not compared; bond lengths are, because rotations
    and relabelling preserve them. `allowed_moved` residue keys are exempt. Returns the largest anchor move."""
    import itertools
    exempt = allowed_moved or set()
    a, b = _atoms(before), _atoms(after)
    problems, worst = [], 0.0
    lost = set(residue_set(before)) - set(residue_set(after))
    if lost:
        problems.append("residues lost: " + ", ".join(f"{c}:{n}" for c, n in sorted(lost)))
    by_res_a, by_res_b = {}, {}
    for k, v in a.items():
        if v[4] not in ("H", "D") and not k[2].startswith(("H", "D")):
            by_res_a.setdefault((k[0], k[1]), {})[k[2]] = v
    for k, v in b.items():
        if v[4] not in ("H", "D") and not k[2].startswith(("H", "D")):
            by_res_b.setdefault((k[0], k[1]), {})[k[2]] = v
    for rk, atoms_a in by_res_a.items():
        if rk in exempt or rk not in by_res_b:
            continue
        atoms_b = by_res_b[rk]
        res = next(iter(atoms_a.values()))[3]
        for name in ANCHORS & atoms_a.keys() & atoms_b.keys():
            d = math.dist(atoms_a[name][:3], atoms_b[name][:3])
            worst = max(worst, d)
            if d > limit:
                problems.append(f"{rk[0]}:{res}{rk[1]}.{name} moved {d:.2f} A")
        for n1, n2 in itertools.combinations(sorted(atoms_a.keys() & atoms_b.keys()), 2):
            if "OXT" in (n1, n2):          # PDB2PQR places the terminal oxygen itself; its bond is not an input bond
                continue
            d0 = math.dist(atoms_a[n1][:3], atoms_a[n2][:3])
            if d0 < BOND_LIMIT:
                d1 = math.dist(atoms_b[n1][:3], atoms_b[n2][:3])
                if abs(d1 - d0) > BOND_CHANGE:
                    problems.append(f"{rk[0]}:{res}{rk[1]} bond {n1}-{n2} changed {d0:.2f} -> {d1:.2f} A")
    if problems:
        _fail(stage, problems)
    return worst


def residues_without_hydrogens(text: str) -> list[tuple]:
    """(chain, number, name) of amino-acid-named residues that carry no hydrogen after protonation: PDB2PQR keeps
    residues it has no template for (MLY, SEC, ASX, UNK, ...) but does not protonate them."""
    per = {}
    for (c, n, name), (x, y, z, res, el) in _atoms(text).items():
        d = per.setdefault((c, n, res), 0)
        per[(c, n, res)] = d + (el in ("H", "D"))
    out = []
    for (c, n, res), h in per.items():
        info = gemmi.find_tabulated_residue(res)
        if h == 0 and (info is None or info.is_amino_acid()) and res not in ("HOH", "WAT"):
            if info is not None or res == "UNK":
                out.append((c, n, res))
    return out


def drop_residues(text: str, keys: list[tuple]) -> str:
    """Remove the given (chain, number, name) residues from PDB text."""
    gone = {(c, n) for c, n, _ in keys}
    return "\n".join(l for l in text.splitlines()
                     if not (l.startswith(("ATOM", "HETATM")) and (l[21], l[22:27].strip()) in gone)) + "\n"


def check_hydrogens(stage: str, text: str):
    """Every amino-acid residue of a protonated model must carry hydrogens, and none may be non-finite."""
    problems = []
    per = {}
    for (c, n, name), (x, y, z, res, el) in _atoms(text).items():
        if not all(map(math.isfinite, (x, y, z))):
            problems.append(f"{c}:{res}{n}.{name} has non-finite coordinates")
        per.setdefault((c, n, res), [0, 0])[0 if el in ("H", "D") else 1] += 1
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
        if el not in ("H", "D"):
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


def missing_amide_hydrogens(text: str) -> list[str]:
    """Residues whose backbone nitrogen carries no hydrogen at all (not Pro): typically the first residue after a
    chain gap, where PDB2PQR places neither the amide H nor the N-terminal ammonium hydrogens. The donor is missing."""
    per: dict = {}
    for (c, n, name), (x, y, z, res, el) in _atoms(text).items():
        d = per.setdefault((c, n), {"res": res, "names": set()})
        d["names"].add(name)
    out = []
    for (c, n), d in per.items():
        info = gemmi.find_tabulated_residue(d["res"])
        if info is None or not info.is_amino_acid() or d["res"] == "PRO" or "N" not in d["names"]:
            continue
        if not d["names"] & {"H", "H1", "H2", "H3", "D", "D1", "D2", "D3"}:
            out.append(f"{c}:{d['res']}{n}")
    return out


# Polar side-chain hydrogens a protonated residue must carry, by residue name (state-specific names for His/Asp/Glu/Lys).
POLAR_SIDECHAIN_H = {
    "ARG": {"HE", "HH11", "HH12", "HH21", "HH22"}, "LYS": {"HZ1", "HZ2", "HZ3"}, "LYN": {"HZ2", "HZ3"},
    "ASN": {"HD21", "HD22"}, "GLN": {"HE21", "HE22"}, "TRP": {"HE1"}, "HID": {"HD1"}, "HIE": {"HE2"},
    "HIP": {"HD1", "HE2"}, "ASH": {"HD2"}, "GLH": {"HE2"}, "SER": {"HG"}, "THR": {"HG1"}, "TYR": {"HH"},
    "CYS": {"HG"},
}


def missing_sidechain_hydrogens(text: str) -> list[str]:
    """'A:ARG134 (HH11)': residues whose heavy atoms are all present but which lack polar hydrogens the residue
    state requires - a donor that is silently absent from the docking receptor."""
    per: dict = {}
    for (c, n, name), (x, y, z, res, el) in _atoms(text).items():
        per.setdefault((c, n), {"res": res, "names": set()})["names"].add(name)
    out = []
    for (c, n), d in per.items():
        need = POLAR_SIDECHAIN_H.get(d["res"])
        if not need:
            continue
        miss = sorted(h for h in need if h not in d["names"])
        if miss:
            out.append(f"{c}:{d['res']}{n} ({','.join(miss)})")
    return out
