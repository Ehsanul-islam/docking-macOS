"""Structure assessment: decide, before any expensive step, whether an input can be docked into.

Every structure gets one verdict:

  ready     nothing to fix
  fixable   the pipeline will repair or exclude things, each named in `findings` (severity "fix")
  unusable  docking into this structure would be meaningless or the pipeline cannot process it
            (severity "block"), with the reason

`findings` also carries "note" items (information the user should see, no action needed). Nothing here
changes the structure; the checks are cheap (a few hundred ms for a typical entry) and deterministic.
"""
from __future__ import annotations

import json
import math
from functools import lru_cache
from pathlib import Path

import gemmi

from . import hetatm
from .errors import InvalidInput

# Backbone bond windows (A). Far outside real chemistry, so only genuinely distorted or mislabelled
# coordinates trigger them; 1960s-era entries still pass unless badly broken.
BACKBONE = {("N", "CA"): (1.25, 1.70), ("CA", "C"): (1.30, 1.75), ("C", "O"): (1.05, 1.45)}
PEPTIDE_BOND = (1.10, 1.75)             # C(i)-N(i+1) when the residues are consecutive and not a gap
MIN_ATOM_DISTANCE = 0.85               # A: no two atoms of one residue may be closer than this
DISTORTED_FRACTION_BLOCK = 0.08         # more than this share of residues distorted => unusable
UNKNOWN_FRACTION_BLOCK = 0.10
TRUNCATED_FRACTION_BLOCK = 0.60         # share of residues missing most side-chain atoms => unusable
CANONICAL = {"ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE", "LEU", "LYS", "MET", "PHE",
             "PRO", "SER", "THR", "TRP", "TYR", "VAL"}
PDB_ATOM_LIMIT = 99_999
LARGE_ATOMS = 150_000
LOW_RESOLUTION = 3.5


@lru_cache(maxsize=1)
def _expected_heavy() -> dict[str, frozenset]:
    import meeko
    t = json.loads((Path(meeko.__file__).parent / "data" / "residue_chem_templates.json").read_text())
    out = {}
    for name, v in t["residue_templates"].items():
        if gemmi.find_tabulated_residue(name) is not None and gemmi.find_tabulated_residue(name).is_amino_acid():
            out[name] = frozenset(a for a in v.get("atom_name", []) if not a.startswith("H"))
    out.setdefault("HIS", out.get("HID", frozenset()))        # Meeko names His by tautomer only
    return out


@lru_cache(maxsize=1)
def _pdbfixer_substitutions() -> frozenset:
    """Modified residues PDBFixer replaces by their canonical parent (MSE -> MET, ...)."""
    try:
        from pdbfixer.pdbfixer import substitutions
        return frozenset(substitutions)
    except Exception:  # noqa: BLE001
        return frozenset(hetatm.MODIFIED)


def truncated_segment_start_prolines(chain) -> list:
    """Prolines that start a chain segment (chain start or after a gap) without their ring atoms CG/CD. PDB2PQR crashes
    on these ("'NoneType' object has no attribute 'coords'": 1UOU, 9NPX); the ring cannot be rebuilt reliably."""
    out, prev = [], None
    for r in chain:
        if not (_is_aa(r) and r.het_flag != "H"):
            prev = None
            continue
        if r.name == "PRO" and {"CG", "CD"} - {a.name for a in r}:
            c, n = (prev.find_atom("C", "*") if prev else None), r.find_atom("N", "*")
            contiguous = prev is not None and not (c and n and c.pos.dist(n.pos) > 2.0)
            if not contiguous:
                out.append(r)
        prev = r
    return out


def _finding(findings, severity, code, message, residues=None):
    f = {"severity": severity, "code": code, "message": message}
    if residues:
        f["residues"] = residues[:50]
        f["n_residues"] = len(residues)
    findings.append(f)


def _is_aa(res) -> bool:
    info = gemmi.find_tabulated_residue(res.name)
    return (info is not None and info.is_amino_acid()) or res.name in hetatm.MODIFIED


def _is_na(res) -> bool:
    info = gemmi.find_tabulated_residue(res.name)
    return info is not None and info.is_nucleic_acid()


def _label(ch, res) -> str:
    return f"{ch.name}:{res.name}{res.seqid.num}{res.seqid.icode.strip()}"


def _distorted(res) -> str | None:
    """Why this residue's own coordinates are impossible, or None."""
    d = {a.name: a for a in res if a.altloc in ("\0", "A", " ")}
    for (a, b), (lo, hi) in BACKBONE.items():
        if a in d and b in d:
            x = d[a].pos.dist(d[b].pos)
            if not lo <= x <= hi:
                return f"{a}-{b} is {x:.2f} A (expected {lo:.2f}-{hi:.2f})"
    # Heavy atoms only: refinement programs place hydroxyl hydrogens 0.84 A from the oxygen, which is normal.
    atoms = [a for a in d.values() if not a.element.is_hydrogen]
    for i, a in enumerate(atoms):
        for b in atoms[i + 1:]:
            x = a.pos.dist(b.pos)
            if x < MIN_ATOM_DISTANCE:
                return f"{a.name} and {b.name} are {x:.2f} A apart"
    return None


def assess(path: str | Path, model: int = 1) -> dict:
    """Assess a receptor file (PDB or mmCIF). Raises InvalidInput only for files that cannot be read at all;
    anything else, including garbage content, is a verdict."""
    from .formats import read_structure
    path = Path(path)
    findings: list[dict] = []
    try:
        st, _ = read_structure(path)
    except InvalidInput as e:
        return _result("unusable", [{"severity": "block", "code": "unreadable", "message": str(e)}], {})
    except Exception as e:  # noqa: BLE001  gemmi raises plain RuntimeError on malformed input
        return _result("unusable", [{"severity": "block", "code": "unreadable",
                                     "message": f"the file could not be parsed: {e}"}], {})
    return assess_structure(st, model, path.suffix.lower())


def assess_structure(st: gemmi.Structure, model: int = 1, suffix: str = "") -> dict:
    """Same as assess() for a structure that is already in memory (e.g. after the user's edits)."""
    findings: list[dict] = []
    if len(st) == 0 or sum(len(ch) for ch in st[0]) == 0:
        return _result("unusable", [{"severity": "block", "code": "empty", "message": "no atoms in the file"}], {})
    if not 1 <= model <= len(st):
        raise InvalidInput(f"model {model} does not exist (the file has {len(st)})")
    m = st[model - 1]
    stats = {"models": len(st), "chains": len(m), "resolution": st.resolution or None}

    method = ""
    try:
        method = (st.info.get("_exptl.method") or "").strip('"\' ')
    except Exception:  # noqa: BLE001
        pass
    stats["method"] = method
    if len(st) > 1:
        _finding(findings, "note", "multiple_models",
                 f"{len(st)} models (NMR ensemble or similar); model {model} is used")
    if st.resolution and st.resolution > LOW_RESOLUTION:
        _finding(findings, "note", "low_resolution",
                 f"resolution {st.resolution:.1f} A: side-chain positions are uncertain")

    n_atoms = sum(len(r) for ch in m for r in ch)
    stats["atoms"] = n_atoms
    if n_atoms > LARGE_ATOMS:
        _finding(findings, "note", "large_structure", f"{n_atoms} atoms: preparation will be slow; "
                 "consider keeping only the chains you need (Delete chain ...)")

    prot, na, unknown = [], [], []
    for ch in m:
        for r in ch:
            if _is_aa(r) and r.het_flag != "H" or (r.name in hetatm.MODIFIED):
                prot.append((ch, r))
            elif _is_na(r):
                na.append((ch, r))
            elif r.name == "UNK":
                unknown.append((ch, r))
    stats.update(protein_residues=len(prot), nucleic_residues=len(na))

    if not prot:
        if na:
            _finding(findings, "block", "nucleic_acid_only",
                     f"no protein: the structure has {len(na)} nucleic-acid residues. Protein-ligand docking "
                     "needs a protein receptor")
        else:
            _finding(findings, "block", "no_protein", "no protein residues found")
        return _result("unusable", findings, stats)
    if na:
        _finding(findings, "note", "nucleic_acid_present",
                 f"{len(na)} nucleic-acid residues present (DNA/RNA is not part of the protein receptor "
                 "unless you choose to keep it)")
    if len(prot) < 20:
        _finding(findings, "note", "short_chain", f"only {len(prot)} protein residues: this is a peptide, not a "
                 "folded protein; results will be unreliable")

    # --- completeness -----------------------------------------------------------------------------
    ca_only = [(c, r) for c, r in prot if sum(1 for a in r if not a.element.is_hydrogen) <= 1]
    if len(ca_only) / len(prot) > 0.9:
        _finding(findings, "block", "ca_only", "Calpha-only model: there are no side chains or backbone atoms "
                 "to dock into")
        return _result("unusable", findings, stats)
    with_side = [(c, r) for c, r in prot if r.name not in ("GLY", "ALA")]
    no_side = [(c, r) for c, r in with_side if r.find_atom("CB", "*") is None]
    if with_side and len(no_side) / len(with_side) > 0.9:
        _finding(findings, "block", "backbone_only", "backbone-only model: side chains are absent from "
                 "almost every residue")
        return _result("unusable", findings, stats)
    unk = [(c, r) for c, r in prot + unknown if r.name == "UNK"]
    if len(unk) / max(1, len(prot) + len(unknown)) > UNKNOWN_FRACTION_BLOCK:
        _finding(findings, "block", "unknown_residues", f"{len(unk)} residues are unidentified (UNK): the "
                 "sequence is unknown, so they cannot be typed")
        return _result("unusable", findings, stats)
    if unk:
        _finding(findings, "fix", "unknown_residues", f"{len(unk)} unidentified (UNK) residues will be left "
                 "out of the receptor", [_label(c, r) for c, r in unk])

    mapped = _pdbfixer_substitutions()
    odd = [(c, r) for c, r in prot if r.name != "UNK" and r.name not in CANONICAL and r.name not in mapped]
    stats["nonstandard_residues"] = len(odd)
    if odd:
        names = sorted({r.name for _, r in odd})
        _finding(findings, "fix", "nonstandard_residues", f"{len(odd)} residues of types without a protonation template "
                 f"({', '.join(names[:6])}) will be left out of the receptor", [_label(c, r) for c, r in odd])
    stubs = [(c, r) for c, r in prot if "CA" not in {a.name for a in r}
             or len({a.name for a in r} & {"N", "C", "O"}) < 2]
    stubs += [(ch, r) for ch in m for r in truncated_segment_start_prolines(ch) if all(r is not x for _, x in stubs)]
    stats["backbone_stubs"] = len(stubs)
    if stubs:
        _finding(findings, "fix", "backbone_stubs", f"{len(stubs)} residues have (almost) no backbone and cannot be "
                 "rebuilt; they will be left out of the receptor", [_label(c, r) for c, r in stubs])
    stub_ids = {id(r) for _, r in stubs}

    expected = _expected_heavy()
    truncated, missing_total = [], 0
    for c, r in prot:
        exp = expected.get(r.name)
        if not exp or id(r) in stub_ids:
            continue
        have = {a.name for a in r}
        miss = exp - have - {"OXT"}
        if miss:
            missing_total += len(miss)
            truncated.append((c, r, len(miss)))
    stats["residues_with_missing_atoms"] = len(truncated)
    if truncated:
        if len(truncated) / len(prot) > TRUNCATED_FRACTION_BLOCK:
            _finding(findings, "block", "mostly_truncated", f"{len(truncated)} of {len(prot)} residues are "
                     "missing side-chain atoms: too much would have to be invented")
            return _result("unusable", findings, stats)
        near_pdb = [f"{_label(c, r)} (-{n})" for c, r, n in truncated]
        _finding(findings, "fix", "missing_atoms", f"{len(truncated)} residues have missing atoms "
                 f"({missing_total} in total); they will be rebuilt by PDB2PQR", near_pdb)

    # --- geometry ---------------------------------------------------------------------------------
    bad = []
    for c, r in prot:
        why = _distorted(r)
        if why:
            bad.append((c, r, why))
    stats["distorted_residues"] = len(bad)
    stats["distorted_keys"] = [f"{c.name}:{r.seqid.num}{r.seqid.icode.strip()}" for c, r, _ in bad]
    if bad:
        share = len(bad) / len(prot)
        lab = [f"{_label(c, r)}: {w}" for c, r, w in bad]
        if share > DISTORTED_FRACTION_BLOCK:
            _finding(findings, "block", "distorted_coordinates",
                     f"distorted coordinates in {len(bad)} of {len(prot)} residues ({100 * share:.0f}%): "
                     "typical of very old or unrefined entries; the result would not be trustworthy", lab)
            return _result("unusable", findings, stats)
        _finding(findings, "fix", "distorted_residues", f"{len(bad)} residues have impossible geometry and "
                 "will be left out of the receptor", lab)

    breaks = []
    for ch in m:
        prev = None
        for r in ch:
            if not (_is_aa(r) and r.het_flag != "H"):
                prev = None
                continue
            if prev is not None:
                c_, n_ = prev.find_atom("C", "*"), r.find_atom("N", "*")
                if c_ and n_ and c_.pos.dist(n_.pos) > 2.0:
                    breaks.append(f"{_label(ch, prev)} -> {_label(ch, r)}")
            prev = r
    stats["chain_breaks"] = len(breaks)
    if breaks:
        _finding(findings, "note", "chain_breaks", f"{len(breaks)} gaps in the chain (unresolved loops); they "
                 "are not modelled unless you ask, and gaps near the site make the pocket incomplete", breaks)

    # --- things the pipeline will do ----------------------------------------------------------------
    alt = [(ch, r) for ch in m for r in ch if any(a.altloc not in ("\0", " ") for a in r)]
    if alt:
        _finding(findings, "fix", "alternate_locations", f"{len(alt)} residues have alternate conformations; "
                 "the highest-occupancy one is kept", [_label(c, r) for c, r in alt])
    n_h = sum(1 for ch in m for r in ch for a in r if a.element.is_hydrogen)
    if n_h:
        _finding(findings, "note", "has_hydrogens", f"{n_h} hydrogens in the input are removed and re-added "
                 "consistently")
    if n_atoms > PDB_ATOM_LIMIT and suffix in (".pdb", ".ent"):
        _finding(findings, "fix", "pdb_atom_limit", f"{n_atoms} atoms exceed the PDB format limit; mmCIF is used "
                 "internally")

    het = {}
    for ch in m:
        for r in ch:
            if r.name in hetatm.WATER or _is_aa(r) and r.het_flag != "H" or _is_na(r):
                continue
            het.setdefault(r.name, [0, sum(1 for a in r if not a.element.is_hydrogen)])[0] += 1
    kinds = {}
    for name, (n, heavy_n) in het.items():
        kinds.setdefault(hetatm.classify(name, heavy_n), []).append(f"{name} x{n}")
    stats["hetero"] = {k: v for k, v in kinds.items()}
    if kinds.get("ligand") or kinds.get("cofactor"):
        _finding(findings, "note", "ligands_present", "ligands/cofactors: " +
                 ", ".join(kinds.get("ligand", []) + kinds.get("cofactor", [])))
    if kinds.get("metal"):
        _finding(findings, "note", "metals_present", "metal ions: " + ", ".join(kinds["metal"]) +
                 " (kept only when within the site radius)")
    covalent = [n for n in het if n in hetatm.GLYCANS]
    if covalent:
        _finding(findings, "note", "glycans_present", "glycans: " + ", ".join(sorted(covalent)))

    verdict = "fixable" if any(f["severity"] == "fix" for f in findings) else "ready"
    return _result(verdict, findings, stats)


def _result(verdict, findings, stats) -> dict:
    reasons = [f["message"] for f in findings if f["severity"] == ("block" if verdict == "unusable" else "fix")]
    return {"verdict": verdict, "summary": "; ".join(reasons) if reasons else "no problems found",
            "findings": findings, "stats": stats}
