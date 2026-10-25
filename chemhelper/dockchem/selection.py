"""Free selection and deletion on a structure (plan §5 stages 3, 5, 6; "Chimera-style" editing).

Selection language (case-insensitive keywords):

  all | none
  chain A,B          resn HOH,GOL,SO*     resi 10-20,45,52A     name CA,CB     elem C,N     altloc A
  protein | nucleic | water | hetero | ligand | cofactor | metal | additive | glycan | modified
  backbone | sidechain | hydrogen        bfactor > 60     occ < 1
  within 5 of <sel>      atoms within 5 A of <sel> (includes <sel> itself)
  around 5 of <sel>      the same, excluding <sel>
  byres <sel>            whole residues that contain any selected atom
  not <sel> | <sel> and <sel> | <sel> or <sel> | ( <sel> )

Editing is non-destructive: the *list of edit commands* is the whole state. Undo is dropping the last command,
redo is putting it back, and replaying the list always gives the same structure. Nothing here mutates the input.
"""
from __future__ import annotations

import fnmatch
import re
from dataclasses import dataclass, field

import gemmi
import numpy as np

from . import hetatm
from .errors import InvalidInput, SelectionError

VARIANTS = {"HID", "HIE", "HIP", "ASH", "GLH", "LYN", "CYX", "CYM"}      # protonation-state names we assign
BACKBONE = {"N", "CA", "C", "O", "OXT", "H", "HA", "H1", "H2", "H3"}
CLASS_WORDS = {"water", "ligand", "cofactor", "metal", "additive", "glycan", "modified"}
KEYWORDS = {"all", "none", "chain", "resn", "resi", "name", "elem", "altloc", "protein", "nucleic", "water",
            "hetero", "ligand", "cofactor", "metal", "additive", "glycan", "modified", "backbone", "sidechain",
            "hydrogen", "bfactor", "occ", "within", "around", "byres", "not", "and", "or", "of"}
PREVIEW_CAP = 500


# ------------------------------------------------------------------------------------------ atom table
@dataclass
class Table:
    n: int
    chain: np.ndarray
    resname: np.ndarray
    seq: np.ndarray
    icode: np.ndarray
    name: np.ndarray
    element: np.ndarray
    altloc: np.ndarray
    b: np.ndarray
    occ: np.ndarray
    xyz: np.ndarray
    het: np.ndarray
    res_id: np.ndarray            # residue group index per atom
    res_class: list[str]          # per residue group: protein|nucleic|water|ligand|...|other
    res_key: list[tuple]          # per residue group: (chain, resname, seq, icode)
    res_slices: list[tuple]       # per residue group: (start, stop) atom range

    @property
    def n_res(self) -> int:
        return len(self.res_key)

    def rid(self, g: int) -> str:
        c, n, s, i = self.res_key[g]
        return f"{c}:{n}{s}{i}"


def _residue_class(name: str, het: bool, n_heavy: int) -> str:
    info = gemmi.find_tabulated_residue(name)
    if name in hetatm.MODIFIED or name in VARIANTS:
        return "protein" if name in VARIANTS else "modified"
    if info is not None and info.is_amino_acid() and not het:
        return "protein"
    if info is not None and info.is_nucleic_acid() and not het:
        return "nucleic"
    c = hetatm.classify(name, n_heavy)
    return {"small_unknown": "ligand"}.get(c, c)


def build_table(st: gemmi.Structure, model: int = 1) -> Table:
    if not 1 <= model <= len(st):
        raise InvalidInput(f"model {model} not in structure (has {len(st)})")
    cols = {k: [] for k in ("chain", "resname", "seq", "icode", "name", "element", "altloc", "b", "occ", "het", "res_id")}
    xyz, res_class, res_key, res_slices = [], [], [], []
    g = 0
    for ch in st[model - 1]:
        for res in ch:
            start = len(cols["name"])
            het = res.het_flag == "H"
            heavy = 0
            for a in res:
                cols["chain"].append(ch.name); cols["resname"].append(res.name); cols["seq"].append(res.seqid.num)
                cols["icode"].append(res.seqid.icode.strip()); cols["name"].append(a.name)
                cols["element"].append(a.element.name.upper()); cols["altloc"].append("" if a.altloc == "\0" else a.altloc)
                cols["b"].append(a.b_iso); cols["occ"].append(a.occ); cols["het"].append(het); cols["res_id"].append(g)
                xyz.append((a.pos.x, a.pos.y, a.pos.z))
                heavy += a.element.name != "H"
            res_class.append(_residue_class(res.name, het, heavy))
            res_key.append((ch.name, res.name, res.seqid.num, res.seqid.icode.strip()))
            res_slices.append((start, len(cols["name"])))
            g += 1
    arr = lambda k, dt=object: np.array(cols[k], dtype=dt)
    return Table(len(cols["name"]), arr("chain"), arr("resname"), arr("seq", np.int64), arr("icode"), arr("name"),
                 arr("element"), arr("altloc"), arr("b", float), arr("occ", float),
                 np.array(xyz, dtype=float).reshape(-1, 3), arr("het", bool), arr("res_id", np.int64),
                 res_class, res_key, res_slices)


# ------------------------------------------------------------------------------------------------ parser
_TOKEN = re.compile(r"\s*(?:(\()|(\))|([<>]=?)|([^\s()<>]+))")


def _tokens(expr: str):
    out, pos = [], 0
    while pos < len(expr):
        m = _TOKEN.match(expr, pos)
        if not m or m.end() == pos:
            break
        tok = m.group(1) or m.group(2) or m.group(3) or m.group(4)
        out.append((tok, m.start(m.lastindex)))
        pos = m.end()
    return out


class _Parser:
    def __init__(self, expr: str):
        self.expr, self.toks, self.i = expr, _tokens(expr), 0

    def err(self, msg, at=None):
        tok_pos = self.toks[self.i][1] if self.i < len(self.toks) else len(self.expr)
        raise SelectionError(msg, tok_pos if at is None else at, self.expr)

    def peek(self):
        return self.toks[self.i][0].lower() if self.i < len(self.toks) else None

    def take(self):
        t = self.toks[self.i][0]
        self.i += 1
        return t

    def parse(self):
        if not self.toks:
            self.err("empty selection")
        node = self.p_or()
        if self.i < len(self.toks):
            self.err(f"unexpected '{self.toks[self.i][0]}'")
        return node

    def p_or(self):
        node = self.p_and()
        while self.peek() == "or":
            self.take(); node = ("or", node, self.p_and())
        return node

    def p_and(self):
        node = self.p_not()
        while self.peek() == "and":
            self.take(); node = ("and", node, self.p_not())
        return node

    def p_not(self):
        if self.peek() == "not":
            self.take(); return ("not", self.p_not())
        return self.p_primary()

    def value(self, what):
        if self.i >= len(self.toks) or self.peek() in ("(", ")", "and", "or", "not", "of"):
            self.err(f"'{what}' needs a value")
        return self.take()

    def p_primary(self):
        t = self.peek()
        if t is None:
            self.err("selection ends unexpectedly")
        if t == "(":
            self.take(); node = self.p_or()
            if self.peek() != ")":
                self.err("missing ')'")
            self.take(); return node
        if t in {"all", "none", "protein", "nucleic", "hetero", "backbone", "sidechain", "hydrogen"} | CLASS_WORDS:
            self.take(); return ("word", t)
        if t in ("chain", "resn", "resi", "name", "elem", "altloc"):
            self.take()
            at = self.toks[self.i][1] if self.i < len(self.toks) else len(self.expr)
            raw = self.value(t)
            vals = [v for v in raw.split(",") if v != ""] or self.err(f"'{t}' needs a value", at)
            if t == "resi":
                for v in vals:                       # validate here so errors are located and need no structure
                    m = re.fullmatch(r"(-?\d+)([A-Za-z]?)(?:-(-?\d+)([A-Za-z]?))?", v)
                    if not m:
                        self.err(f"bad residue range '{v}' (use 10, 10-20 or 52A)", at)
                    if m.group(3) is not None and int(m.group(3)) < int(m.group(1)):
                        self.err(f"range '{v}' runs backwards", at)
            return (t, vals)
        if t in ("bfactor", "occ"):
            self.take()
            op = self.peek()
            if op not in ("<", ">", "<=", ">="):
                self.err(f"'{t}' needs a comparison, e.g. {t} > 50")
            self.take()
            try:
                return (t, op, float(self.take()))
            except (ValueError, IndexError):
                self.err(f"'{t}' needs a number")
        if t in ("within", "around"):
            self.take()
            try:
                dist = float(self.take())
            except (ValueError, IndexError):
                self.err(f"'{t}' needs a distance in angstrom, e.g. {t} 5 of ligand")
            if dist <= 0 or dist > 50:
                self.err("distance must be between 0 and 50 angstrom")
            if self.peek() != "of":
                self.err(f"expected 'of' after the distance in '{t}'")
            self.take()
            return (t, dist, self.p_not())
        if t == "byres":
            self.take(); return ("byres", self.p_not())
        self.err(f"unknown keyword '{self.toks[self.i][0]}'")


def parse(expr: str):
    return _Parser(expr).parse()


# --------------------------------------------------------------------------------------------- evaluator
def _resi_match(t: Table, specs: list[str], pos: int) -> np.ndarray:
    mask = np.zeros(t.n, bool)
    for s in specs:
        m = re.fullmatch(r"(-?\d+)([A-Za-z]?)(?:-(-?\d+)([A-Za-z]?))?", s)
        if not m:
            raise SelectionError(f"bad residue range '{s}' (use 10, 10-20 or 52A)", pos)
        lo, ic, hi, _ = m.groups()
        if hi is None:
            mask |= (t.seq == int(lo)) & ((t.icode == ic.upper()) if ic else True)
        else:
            if int(hi) < int(lo):
                raise SelectionError(f"range '{s}' runs backwards", pos)
            mask |= (t.seq >= int(lo)) & (t.seq <= int(hi))
    return mask


def _glob(arr: np.ndarray, patterns: list[str], upper=True) -> np.ndarray:
    pats = [p.upper() if upper else p for p in patterns]
    exact = {p for p in pats if not any(c in p for c in "*?[")}
    wild = [p for p in pats if p not in exact]
    mask = np.isin(arr, list(exact)) if exact else np.zeros(len(arr), bool)
    for w in wild:
        mask |= np.array([fnmatch.fnmatchcase(x, w) for x in arr], bool) if len(arr) else mask
    return mask


def evaluate(t: Table, node, expr: str = "") -> np.ndarray:
    kind = node[0]
    if kind == "word":
        w = node[1]
        if w == "all":
            return np.ones(t.n, bool)
        if w == "none":
            return np.zeros(t.n, bool)
        if w == "hetero":
            return t.het.copy()
        if w == "hydrogen":
            return t.element == "H"
        if w == "backbone":
            return np.isin(t.name, list(BACKBONE)) & np.array([t.res_class[g] == "protein" for g in t.res_id], bool)
        if w == "sidechain":
            return ~np.isin(t.name, list(BACKBONE)) & (t.element != "H") & \
                np.array([t.res_class[g] == "protein" for g in t.res_id], bool)
        cls = np.array([c == w for c in t.res_class], bool)
        return cls[t.res_id] if t.n else np.zeros(0, bool)
    if kind == "chain":
        return np.isin(t.chain, node[1])
    if kind == "resn":
        return _glob(t.resname, node[1])
    if kind == "name":
        return _glob(t.name, node[1])
    if kind == "elem":
        return np.isin(t.element, [e.upper() for e in node[1]])
    if kind == "altloc":
        return np.isin(t.altloc, [a if a != "." else "" for a in node[1]])
    if kind == "resi":
        return _resi_match(t, node[1], 0)
    if kind in ("bfactor", "occ"):
        arr = t.b if kind == "bfactor" else t.occ
        op, v = node[1], node[2]
        return {"<": arr < v, ">": arr > v, "<=": arr <= v, ">=": arr >= v}[op]
    if kind == "not":
        return ~evaluate(t, node[1], expr)
    if kind == "and":
        return evaluate(t, node[1], expr) & evaluate(t, node[2], expr)
    if kind == "or":
        return evaluate(t, node[1], expr) | evaluate(t, node[2], expr)
    if kind == "byres":
        sel = evaluate(t, node[1], expr)
        hit = np.zeros(t.n_res, bool)
        hit[np.unique(t.res_id[sel])] = True
        return hit[t.res_id] if t.n else sel
    if kind in ("within", "around"):
        from scipy.spatial import cKDTree
        target = evaluate(t, node[2], expr)
        if not target.any() or not t.n:
            return np.zeros(t.n, bool)
        near = cKDTree(t.xyz).query_ball_point(t.xyz[target], r=node[1])
        mask = np.zeros(t.n, bool)
        idx = np.unique(np.concatenate([np.asarray(x, dtype=np.int64) for x in near])) if len(near) else []
        mask[idx] = True
        return mask & ~target if kind == "around" else mask
    raise SelectionError(f"cannot evaluate '{kind}'", 0, expr)


def select(t: Table, expr: str) -> np.ndarray:
    return evaluate(t, parse(expr), expr)


# ------------------------------------------------------------------------------------------ edit session
@dataclass
class Edit:
    command: str                    # "delete <sel>" | "keep <sel>"
    n_atoms: int = 0
    n_residues: int = 0
    residues: list[str] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)


def _split(command: str) -> tuple[str, str]:
    m = re.match(r"\s*(delete|keep|remove)\s+(.*\S)\s*$", command, re.I | re.S)
    if not m:
        raise SelectionError("an edit must look like 'delete <selection>' or 'keep <selection>'", 0, command)
    return ("delete" if m.group(1).lower() in ("delete", "remove") else "keep"), m.group(2)


def _breaks(t: Table, alive: np.ndarray) -> int:
    """Chain breaks among alive protein residues, by C(i)-N(i+1) distance (works without SEQRES)."""
    n, prev = 0, None
    for g in range(t.n_res):
        s, e = t.res_slices[g]
        if t.res_class[g] != "protein" or not alive[s:e].any():
            prev = None if t.res_class[g] != "protein" else prev
            continue
        names = {t.name[i]: i for i in range(s, e) if alive[i]}
        if prev is not None and prev[0] == t.chain[s]:
            if "C" in prev[1] and "N" in names and np.linalg.norm(t.xyz[prev[1]["C"]] - t.xyz[names["N"]]) > 2.0:
                n += 1
        prev = (t.chain[s], names)
    return n


class EditSession:
    """Applies a list of edit commands to a structure, non-destructively."""

    def __init__(self, st: gemmi.Structure, model: int = 1):
        self.st, self.model = st, model
        self.table = build_table(st, model)
        self.alive = np.ones(self.table.n, bool)
        self.edits: list[Edit] = []

    # -- previewing never changes state
    def preview(self, command: str) -> Edit:
        verb, expr = _split(command)
        t = self.table
        sel = select(t, expr) & self.alive
        removed = sel if verb == "delete" else (~sel & self.alive)
        return self._describe(command, removed)

    def _describe(self, command: str, removed: np.ndarray) -> Edit:
        t = self.table
        groups = np.unique(t.res_id[removed]) if removed.any() else np.array([], int)
        partial = [t.rid(g) for g in groups
                   if removed[t.res_slices[g][0]:t.res_slices[g][1]].sum()
                   < self.alive[t.res_slices[g][0]:t.res_slices[g][1]].sum()]
        warns = []
        if partial:
            warns.append(f"{len(partial)} residue(s) would lose only some atoms (e.g. {partial[0]}); partial residues "
                         "cannot be prepared. Use byres to select whole residues.")
        prot_total = sum(1 for g in range(t.n_res) if t.res_class[g] == "protein" and self.alive[t.res_slices[g][0]:t.res_slices[g][1]].any())
        prot_removed = sum(1 for g in groups if t.res_class[g] == "protein")
        if prot_total and prot_removed == prot_total:
            warns.append("this removes every protein residue")
        elif prot_total and prot_removed / prot_total > 0.5:
            warns.append(f"this removes {prot_removed} of {prot_total} protein residues")
        after = self.alive & ~removed
        new_breaks = _breaks(t, after) - _breaks(t, self.alive)
        if new_breaks > 0:
            warns.append(f"this creates {new_breaks} new chain break(s); residues next to a break may fail to prepare")
        lig = np.array([t.res_class[g] == "ligand" for g in t.res_id], bool) & after
        if lig.any() and removed.any():
            from scipy.spatial import cKDTree
            pocket = cKDTree(t.xyz[lig]).query_ball_point(t.xyz[removed & (t.element != "H")], r=5.0)
            hit = {int(t.res_id[i]) for i, near in zip(np.flatnonzero(removed & (t.element != "H")), pocket) if near}
            hit = [g for g in hit if t.res_class[g] == "protein"]
            if hit:
                warns.append(f"{len(hit)} residue(s) within 5 A of a ligand would be removed (pocket residues)")
        return Edit(command, int(removed.sum()), len(groups), [t.rid(g) for g in groups[:PREVIEW_CAP]], warns)

    # -- applying
    def apply(self, command: str) -> Edit:
        e = self.preview(command)
        verb, expr = _split(command)
        sel = select(self.table, expr) & self.alive
        self.alive = self.alive & ~(sel if verb == "delete" else ~sel)
        self.edits.append(e)
        return e

    def apply_all(self, commands: list[str]) -> list[Edit]:
        return [self.apply(c) for c in commands]

    def commands(self) -> list[str]:
        return [e.command for e in self.edits]

    # -- output
    def summary(self) -> dict:
        t, alive = self.table, self.alive
        groups = np.unique(t.res_id[alive]) if alive.any() else []
        by_class: dict[str, int] = {}
        for g in groups:
            by_class[t.res_class[g]] = by_class.get(t.res_class[g], 0) + 1
        return {"atoms": int(alive.sum()), "residues": len(groups), "by_class": by_class,
                "chains": sorted({str(c) for c in t.chain[alive]})}

    def materialize(self) -> gemmi.Structure:
        """A new structure with only the surviving atoms (other models are dropped)."""
        out = self.st.clone()
        for i in range(len(out) - 1, -1, -1):          # keep only the model the table was built from
            if i != self.model - 1:
                del out[i]
        m = out[0]
        k = len(self.alive)
        assert sum(len(res) for ch in m for res in ch) == k, "table and structure disagree"
        pos = k
        for ci in range(len(m) - 1, -1, -1):
            ch = m[ci]
            for ri in range(len(ch) - 1, -1, -1):
                res = ch[ri]
                for ai in range(len(res) - 1, -1, -1):
                    pos -= 1
                    if not self.alive[pos]:
                        del res[ai]
                if len(res) == 0:
                    del ch[ri]
            if len(ch) == 0:
                del m[ci]
        out.setup_entities()
        return out
