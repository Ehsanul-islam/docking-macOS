"""Per-residue protonation overrides and Asn/Gln/His flip decisions (gap 4).

PDB2PQR + PROPKA decide protonation states and flip Asn/Gln/His side chains to optimise hydrogen
bonds. Both are judgement calls a user may want to overrule (a catalytic Asp that is protonated, a His
whose tautomer is known from the literature, a flip that breaks a contact seen in the crystal).

Overrides are applied by naming the residue (HIE, ASH, ...) in the PDB2PQR input, which PDB2PQR
respects; the result is then *verified* against the hydrogens actually present, so an override that
did not take effect is an error, never a silent no-op. Flips are detected by comparing PDB2PQR's
output with its input; a rejected flip is undone by rotating the side-chain amide / imidazole group
back by 180 degrees about its axis.
"""
from __future__ import annotations

import math

from .errors import InvalidInput, ToolFailure

# Allowed states per parent residue. CYX is deliberately absent: disulfides are detected from the
# SG-SG distance, not chosen by the user.
OVERRIDE_STATES = {
    "HIS": {"HID", "HIE", "HIP"},
    "ASP": {"ASP", "ASH"},
    "GLU": {"GLU", "GLH"},
    "LYS": {"LYS", "LYN"},
    "CYS": {"CYS", "CYM"},
}
_PARENT_OF = {s: p for p, ss in OVERRIDE_STATES.items() for s in ss}

# Hydrogens that must be present / absent for a state to be genuinely what was asked for.
_H_PRESENT = {"HID": {"HD1"}, "HIE": {"HE2"}, "HIP": {"HD1", "HE2"}, "ASH": {"HD2"}, "GLH": {"HE2"},
              "LYN": {"HZ2", "HZ3"}, "LYS": {"HZ1", "HZ2", "HZ3"}}
_H_ABSENT = {"HID": {"HE2"}, "HIE": {"HD1"}, "ASP": {"HD2"}, "GLU": {"HE2"}, "LYN": {"HZ1"}, "CYM": {"HG"}}

# Flip groups: (axis atoms, atoms rotated). Asn/Gln flip swaps O and N; His swaps the ring.
FLIP_GROUPS = {
    "ASN": (("CB", "CG"), ("OD1", "ND2", "HD21", "HD22")),
    "GLN": (("CG", "CD"), ("OE1", "NE2", "HE21", "HE22")),
    "HIS": (("CB", "CG"), ("ND1", "CD2", "CE1", "NE2", "HD1", "HD2", "HE1", "HE2")),
}
_FLIP_PAIRS = {"ASN": (("OD1", "ND2"),), "GLN": (("OE1", "NE2"),), "HIS": (("ND1", "CD2"), ("NE2", "CE1"))}
for _s in ("HID", "HIE", "HIP"):
    FLIP_GROUPS[_s] = FLIP_GROUPS["HIS"]
    _FLIP_PAIRS[_s] = _FLIP_PAIRS["HIS"]
FLIP_AXIS_TOLERANCE = 0.3   # A: the atoms the group rotates about must not have moved
FLIP_MIN_MOVE = 1.0         # A: mean heavy-atom movement of the group that counts as a flip


def parse_key(key: str) -> tuple[str, str]:
    """'A:45' or 'A:45B' -> ('A', '45B'). Residue names are not part of the key: the structure decides."""
    chain, sep, num = str(key).partition(":")
    num = num.strip()
    if not sep or len(chain) != 1 or not num or not num.lstrip("-")[:1].isdigit():
        raise InvalidInput(f"bad residue key {key!r}: use chain:number, e.g. 'A:45' or 'A:45B'")
    return chain, num


def _residue_lines(text: str):
    """{(chain, number): [line indices]} for ATOM records."""
    out: dict = {}
    for i, l in enumerate(text.splitlines()):
        if l.startswith(("ATOM", "HETATM")):
            out.setdefault((l[21], l[22:27].strip()), []).append(i)
    return out


def _name(l: str) -> str:
    return l[17:20].strip()


def apply_overrides(pdb_text: str, overrides: dict[str, str]) -> tuple[str, list[dict]]:
    """Rename residues to the requested protonation state. Returns the new text and one record per
    override ({key, old, new}); raises InvalidInput for unknown residues or impossible states."""
    if not overrides:
        return pdb_text, []
    lines = pdb_text.splitlines()
    index = _residue_lines(pdb_text)
    done = []
    for raw_key, state in overrides.items():
        key = parse_key(raw_key)
        state = str(state).strip().upper()
        if key not in index:
            raise InvalidInput(f"protonation override {raw_key}: no such residue in the receptor "
                               "(it may have been deleted by an edit, or belongs to an unselected chain)")
        old = _name(lines[index[key][0]])
        parent = "HIS" if old == "HIS" else _PARENT_OF.get(old)
        if parent is None:
            raise InvalidInput(f"protonation override {raw_key}: {old} has no selectable protonation "
                               f"state (supported: {', '.join(sorted(OVERRIDE_STATES))})")
        if state not in OVERRIDE_STATES[parent]:
            raise InvalidInput(f"protonation override {raw_key}: {state} is not a state of {parent} "
                               f"(choose from {', '.join(sorted(OVERRIDE_STATES[parent]))})")
        for i in index[key]:
            lines[i] = lines[i][:17] + state.ljust(3) + lines[i][20:]
        done.append({"key": f"{key[0]}:{key[1]}", "old": old, "new": state})
    return "\n".join(lines) + "\n", done


def verify_overrides(prot_text: str, applied: list[dict]) -> None:
    """The override must have survived PDB2PQR with the matching hydrogens, otherwise fail loudly."""
    index = _residue_lines(prot_text)
    lines = prot_text.splitlines()
    for a in applied:
        key = parse_key(a["key"])
        rows = index.get(key)
        if not rows:
            raise ToolFailure(f"protonation override {a['key']}: residue missing from PDB2PQR output")
        got = _name(lines[rows[0]])
        names = {lines[i][12:16].strip() for i in rows}
        want = a["new"]
        if got != want:
            raise ToolFailure(f"protonation override {a['key']} -> {want} was not honoured "
                              f"(PDB2PQR produced {got})")
        missing = _H_PRESENT.get(want, set()) - names
        extra = _H_ABSENT.get(want, set()) & names
        if missing or extra:
            raise ToolFailure(f"protonation override {a['key']} -> {want}: hydrogens do not match the "
                              f"requested state (missing {sorted(missing)}, unexpected {sorted(extra)})")


def _coords(text: str):
    out: dict = {}
    for l in text.splitlines():
        if l.startswith(("ATOM", "HETATM")):
            out.setdefault((l[21], l[22:27].strip()), {})[l[12:16].strip()] = (
                float(l[30:38]), float(l[38:46]), float(l[46:54]), _name(l))
    return out


def _d(a, b) -> float:
    return math.dist(a[:3], b[:3])


def detect_flips(before_text: str, after_text: str) -> list[dict]:
    """Residues whose terminal group PDB2PQR rotated about its axis: the axis atoms stay put while the
    group's heavy atoms move by FLIP_MIN_MOVE on average (the groups are only roughly symmetric, so an
    exact atom-swap test misses real flips)."""
    b, a = _coords(before_text), _coords(after_text)
    flips = []
    for key, atoms_a in a.items():
        atoms_b = b.get(key)
        if not atoms_b:
            continue
        name = next(iter(atoms_a.values()))[3]
        group = FLIP_GROUPS.get(name)
        if not group:
            continue
        axis, moving = group
        heavy = [n for n in moving if not n.startswith("H") and n in atoms_a and n in atoms_b]
        if len(heavy) < 2 or not all(n in atoms_a and n in atoms_b for n in axis):
            continue
        if any(_d(atoms_a[n], atoms_b[n]) > FLIP_AXIS_TOLERANCE for n in axis):
            continue
        mean = sum(_d(atoms_a[n], atoms_b[n]) for n in heavy) / len(heavy)
        if mean >= FLIP_MIN_MOVE:
            flips.append({"key": f"{key[0]}:{key[1]}", "resname": name})
    return flips


def _rotate_half_turn(p, a0, a1):
    ax = [a1[i] - a0[i] for i in range(3)]
    n = math.sqrt(sum(c * c for c in ax))
    ax = [c / n for c in ax]
    v = [p[i] - a0[i] for i in range(3)]
    t = sum(v[i] * ax[i] for i in range(3))
    proj = [a0[i] + t * ax[i] for i in range(3)]
    return [2 * proj[i] - p[i] for i in range(3)]


def undo_flips(pdb_text: str, keys: list[str]) -> str:
    """Rotate the flipped group of each residue back by 180 degrees about its axis (hydrogens included)."""
    if not keys:
        return pdb_text
    lines = pdb_text.splitlines()
    index = _residue_lines(pdb_text)
    for raw in keys:
        key = parse_key(raw)
        rows = index.get(key)
        if not rows:
            raise InvalidInput(f"flip decision {raw}: no such residue")
        name = _name(lines[rows[0]])
        group = FLIP_GROUPS.get(name)
        if not group:
            raise InvalidInput(f"flip decision {raw}: {name} is not an Asn, Gln or His")
        axis, moving = group
        by_name = {lines[i][12:16].strip(): i for i in rows}
        if not all(n in by_name for n in axis):
            raise ToolFailure(f"flip decision {raw}: axis atoms missing")
        pt = lambda i: [float(lines[i][30:38]), float(lines[i][38:46]), float(lines[i][46:54])]  # noqa: E731
        a0, a1 = pt(by_name[axis[0]]), pt(by_name[axis[1]])
        for n in moving:
            if n in by_name:
                i = by_name[n]
                x, y, z = _rotate_half_turn(pt(i), a0, a1)
                lines[i] = lines[i][:30] + f"{x:8.3f}{y:8.3f}{z:8.3f}" + lines[i][54:]
    return "\n".join(lines) + "\n"
