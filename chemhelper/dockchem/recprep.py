"""Receptor preparation pipeline (plan §5).

gemmi (structure edits) -> PDBFixer (missing atoms, non-standard residues)
-> PDB2PQR + PROPKA (protonation, flips) -> [OpenMM restrained minimization] -> Meeko (PDBQT).
Every change is recorded in the Prep Report (§5.2). Kept HETATM groups bypass PDBFixer/PDB2PQR and
are merged back before Meeko (they get no hydrogens: a documented limitation).
"""
from __future__ import annotations

import io
import math
import re
import subprocess
import sys
import tempfile
from dataclasses import dataclass, field
from pathlib import Path

import gemmi

from . import cancel, hetatm, protonation, resources
from . import stagecheck
from .assess import assess_structure, truncated_segment_start_prolines
from .errors import FormatLimit, HelperError, InvalidInput, ToolFailure, UnusableStructure
from .provenance import record
from . import __version__

LOOP_NEAR_SITE = 10.0
CLASH_CUTOFF = 2.2                 # heavy-atom non-bonded contact, angstrom
RESTRAINT_K = 41840.0            # 10 kcal/mol/A^2 in kJ/mol/nm^2 (user-requested minimisation)
REPAIR_K = 418400.0              # 100 kcal/mol/A^2: hydrogen repair must hold heavy atoms nearly fixed
UNUSUAL_STATES = {"ASH": "protonated Asp", "GLH": "protonated Glu", "LYN": "neutral Lys",
                  "CYM": "deprotonated Cys", "HIP": "doubly protonated His"}
TERMINAL_RENAME = {"HID", "HIE", "HIP", "ASH", "GLH", "LYN", "CYX", "CYM"}
CHARGE = {"ASP": -1, "GLU": -1, "LYS": 1, "ARG": 1, "HIP": 1, "CYM": -1}


@dataclass
class ReceptorPrepOptions:
    model: int = 1
    chains: list[str] | None = None            # None = auto
    reference_ligand: str | None = None        # "LIG" or "A:LIG:101"
    site_radius: float = 5.0
    water_keep_radius: float | None = None     # None = remove all waters
    hetatm_overrides: dict = field(default_factory=dict)   # resname -> class
    nonstandard: str = "map"                   # map | remove
    # Who rebuilds missing heavy atoms: PDB2PQR (default), PDBFixer, or nobody. PDBFixer's rebuild of
    # heavily truncated residues was geometrically unreliable on real entries (1CLL, 1EMA).
    rebuild_missing_atoms: str = "pdb2pqr"      # pdb2pqr | pdbfixer | none
    model_loops_max: int = 0                   # model internal gaps up to N residues (0 = flag only)
    ph: float = 7.4
    optimize_flips: bool = True
    # Per-residue protonation, e.g. {"A:64": "HIE", "A:25": "ASH"} (see protonation.py). Applied before
    # PDB2PQR and verified afterwards; an override that does not take effect is an error.
    protonation_overrides: dict = field(default_factory=dict)
    # Per-residue Asn/Gln/His flip decisions, e.g. {"A:45": "reject"}. Every flip PDB2PQR makes is listed
    # in the report (stage "flip"); "reject" restores the crystal orientation. Default is accept.
    flip_decisions: dict = field(default_factory=dict)
    minimize: bool = False
    flex_residues: list[str] = field(default_factory=list)   # e.g. ["A:15"]
    # Free edits applied before anything else, e.g. ["delete chain B", "delete resn SO4,GOL",
    # "delete within 3 of resn XYZ and water"] (see selection.py). Each one is logged in the report.
    edits: list[str] = field(default_factory=list)
    # Who writes the PDBQT. "auto" (default): the built-in table writer (rigid only, see recwriter.py), which
    # cannot be broken by odd geometry; Meeko steps in when the table lacks rules for some residue, and
    # otherwise runs as an independent cross-check (disagreements are reported, never fatal).
    pdbqt_writer: str = "auto"                  # auto | table | meeko
    crosscheck: bool = True
    # Structures the assessment calls unusable are refused with the reasons; True prepares them anyway.
    allow_unusable: bool = False
    alphafold: bool = False
    min_plddt: float = 70.0
    trim_low_plddt_termini: bool = False


@dataclass
class ReceptorResult:
    pdb: Path
    pdbqt: Path | None
    flex_pdbqt: Path | None
    ligands: list[dict]
    report: dict


class _Report:
    def __init__(self):
        self.changes: list[dict] = []
        self.warnings: list[str] = []
        self.tools: list[dict] = []
        self.dropped: list[dict] = []        # residues absent from the PDBQT: makes the result "degraded"
        self.pdbqt_info: dict = {}           # which writer produced the PDBQT and how the cross-check went
        self.pending_het: dict = {}          # (chain, num, name) -> PDB lines Meeko has no template for

    def change(self, stage, residue, change, reason, tool):
        self.changes.append({"stage": stage, "residue": residue, "change": change,
                             "reason": reason, "tool": tool})

    def warn(self, msg):
        if msg not in self.warnings:
            self.warnings.append(msg)


# Messages that record a routine, expected step rather than a possible problem for the docking result. They stay in
# `warnings` (nothing is hidden) and are tagged "notice" in `warning_details` so a UI can fold them away.
NOTICE_PREFIXES = ("chains renamed", "unusual protonation", "metal ", "pocket residue", "flip decision")


def warning_severity(msg: str) -> str:
    return "notice" if msg.startswith(NOTICE_PREFIXES) else "warning"


def _rid(chain, res) -> str:
    return f"{chain}:{res.name}{res.seqid.num}{res.seqid.icode.strip()}"


def _heavy(res) -> int:
    return sum(1 for a in res if not a.element.is_hydrogen)


def _is_polymer(res) -> bool:
    if res.name in hetatm.MODIFIED:
        return True
    info = gemmi.find_tabulated_residue(res.name)
    # A HETATM record with a standard amino-acid name is a bound peptide (e.g. the Gly-Tyr inhibitor in
    # 3CPA), not part of the receptor chain.
    return info is not None and info.is_amino_acid() and res.het_flag != "H"


def _min_dist(res, pos_list) -> float:
    best = math.inf
    for a in res:
        for p in pos_list:
            best = min(best, a.pos.dist(p))
    return best


def _collapse_altlocs(model, rep):
    for ch in model:
        for res in ch:
            alts = []
            for a in res:
                if a.altloc != "\0" and a.altloc not in alts:
                    alts.append(a.altloc)
            if len(alts) < 2:
                for a in res:
                    a.altloc = "\0"
                continue
            occ = {alt: sum(a.occ for a in res if a.altloc == alt) for alt in alts}
            keep = max(alts, key=lambda k: (occ[k], -alts.index(k)))     # first on tie
            for name, alt in [(a.name, a.altloc) for a in res if a.altloc not in ("\0", keep)]:
                res.remove_atom(name, alt)
            for a in res:
                a.altloc = "\0"
            rep.change("altlocs", _rid(ch.name, res), f"kept altloc {keep} of {''.join(alts)}",
                       "highest occupancy (first on tie)", "gemmi")


ALIAS_POOL = [chr(c) for c in list(range(65, 91)) + list(range(97, 123)) + list(range(48, 58))]


def _alias_chains(st, opts, rep):
    """The PDB format (used by PDBFixer, PDB2PQR and Meeko) holds one-letter chain names; mmCIF entries can have
    'AAA' or 'AB1'. Rename such chains to free letters and tell the user the mapping. Returns (opts, map)."""
    import dataclasses
    model = st[0]
    long_ = [ch.name for ch in model if len(ch.name) > 1]
    if not long_:
        return opts, {}
    free = [c for c in ALIAS_POOL if c not in {ch.name for ch in model}]
    if len(long_) > len(free) and opts.chains is not None:
        # more long names than letters: only the chains the user picked need a one-letter name, the rest are left out
        unused = [c for c in long_ if c not in opts.chains]
        for name in unused:
            model.remove_chain(name)
            rep.change("chains", name, "left out", "multi-letter chain name does not fit the PDB format and the "
                       "chain was not selected", "dockchem")
        long_ = [c for c in long_ if c in opts.chains]
    if len(long_) > len(free):
        raise InvalidInput(f"{len(long_)} chains with multi-letter names but only {len(free)} free one-letter names: "
                           "choose the chains to dock with the chains option (e.g. chains=['AAA']) so the rest can be left out")
    mapping = dict(zip(long_, free))
    for ch in model:
        if ch.name in mapping:
            ch.name = mapping[ch.name]
    for old, new in mapping.items():
        rep.change("chains", old, f"renamed to {new}", "multi-letter chain names do not fit the PDB format", "dockchem")
    rep.warn("chains renamed for processing: " + ", ".join(f"{o}->{n}" for o, n in mapping.items()))

    def tr(spec):
        if not spec:
            return spec
        head, sep, rest = spec.partition(":")
        return mapping.get(head, head) + sep + rest if sep else spec
    new_opts = dataclasses.replace(
        opts, chains=[mapping.get(c, c) for c in opts.chains] if opts.chains is not None else None,
        reference_ligand=tr(opts.reference_ligand), flex_residues=[tr(f) for f in opts.flex_residues])
    return new_opts, mapping


def _is_backbone_stub(res) -> str | None:
    """PDB2PQR cannot rebuild a residue that has (almost) no backbone ('Too few atoms present to reconstruct')."""
    have = {a.name for a in res}
    if "CA" not in have or len(have & {"N", "C", "O"}) < 2:
        return "backbone atoms missing (has " + ",".join(sorted(have - {"OXT"}))[:20] + ")"
    return None


AD4_METALS = {"Mg": "Mg", "Ca": "Ca", "Mn": "Mn", "Fe": "Fe", "Zn": "Zn"}   # Meeko's AD4 typing (ad4_types.json)
ORGANIC = {"C", "N", "O", "S", "P", "F", "Cl", "Br", "I", "B"}
SITE_IMPACT_RADIUS = 8.0     # A: something left out this close to the site can change a docking result


def _unsupported_reason(res) -> str | None:
    """Why a kept hetero group cannot be written to the receptor PDBQT, or None if it can."""
    els = {a.element.name for a in res if not a.element.is_hydrogen}
    odd = els - ORGANIC
    if not odd:
        return None
    if len(res) == 1:
        el = next(iter(els))
        return None if el in AD4_METALS else f"{el} has no AutoDock atom type"
    return f"{res.name} is a metal complex ({', '.join(sorted(odd))}); metal cofactors cannot be represented yet"


def _dropped_kind(label: str, reason: str) -> str:
    """glycan | residue_defect (stub, distorted, no template) | cofactor_or_ion (could be part of a pocket)."""
    name = re.sub(r"-?\d+[A-Za-z]?$", "", label.partition(":")[2])
    if name in hetatm.GLYCANS:
        return "glycan"
    if reason.startswith(("backbone atoms missing", "impossible", "no protonation template")):
        return "residue_defect"
    return "cofactor_or_ion"


def _annotate_dropped(rep, model, site_pos) -> None:
    """Add distance_to_site / at_site / kind to every left-out residue. With a site, distance decides. Without one
    (blind docking) nothing can be located: sugars and single defective residues are only noted, but a left-out
    cofactor or ion could be part of a pocket, so it still counts."""
    index = {}
    for ch in model:
        for res in ch:
            index.setdefault((ch.name, f"{res.seqid.num}{res.seqid.icode.strip()}"), []).append(res)
    for d in rep.dropped:
        chain, _, rest = d["residue"].partition(":")
        m = re.search(r"(-?\d+[A-Za-z]?)$", rest)
        found = index.get((chain, m.group(1))) if m else None
        d["kind"] = _dropped_kind(d["residue"], d["reason"])
        if site_pos and found:
            dist = min(_min_dist(r, site_pos) for r in found)
            d["distance_to_site"] = round(dist, 1)
            d["at_site"] = dist <= SITE_IMPACT_RADIUS
        else:
            d["distance_to_site"], d["at_site"] = None, None


def _is_degrading(d: dict) -> bool:
    if d.get("at_site") is True:
        return True
    return d.get("at_site") is None and d.get("kind") == "cofactor_or_ion"


def _find_reference(model, spec):
    if not spec:
        return []
    parts = spec.split(":")
    out = []
    for ch in model:
        for res in ch:
            if res.het_flag != "H":
                continue
            if len(parts) == 1 and res.name == parts[0]:
                out.append((ch.name, res))
            elif len(parts) == 3 and (ch.name, res.name, str(res.seqid.num)) == tuple(parts):
                out.append((ch.name, res))
    return out


def _chain_breaks(model):
    """Geometric chain-break detection (works without SEQRES)."""
    breaks = []
    for ch in model:
        prev = None
        for res in ch:
            if not _is_polymer(res):
                prev = None
                continue
            if prev is not None:
                c, n = prev.find_atom("C", "*"), res.find_atom("N", "*")
                if c and n and c.pos.dist(n.pos) > 2.0:
                    breaks.append((ch.name, prev, res))
            prev = res
    return breaks


def _alphafold(st_model, site_res_ids, opts, rep):
    for ch in st_model:
        low = []
        for res in ch:
            if _is_polymer(res):
                ca = res.find_atom("CA", "*")
                if ca and ca.b_iso < opts.min_plddt:
                    low.append(res)
                    if _rid(ch.name, res) in site_res_ids:
                        rep.warn(f"pocket residue {_rid(ch.name, res)} has pLDDT {ca.b_iso:.0f} < {opts.min_plddt:g}")
        if opts.trim_low_plddt_termini:
            for end in ("start", "end"):
                while len(ch) and True:
                    res = ch[0] if end == "start" else ch[len(ch) - 1]
                    ca = res.find_atom("CA", "*")
                    if _is_polymer(res) and ca and ca.b_iso < opts.min_plddt:
                        rep.change("alphafold", _rid(ch.name, res), f"trimmed {end} residue",
                                   f"pLDDT {ca.b_iso:.0f} < {opts.min_plddt:g}", "gemmi")
                        del ch[0 if end == "start" else len(ch) - 1]
                    else:
                        break


# Helper tools (Meeko, PDB2PQR) can fetch definitions from RCSB on their own. Local-first (plan principle 5)
# means they must never do so silently, so they run with an unroutable proxy.
OFFLINE_ENV = {"HTTP_PROXY": "http://127.0.0.1:9", "HTTPS_PROXY": "http://127.0.0.1:9",
               "ALL_PROXY": "http://127.0.0.1:9", "http_proxy": "http://127.0.0.1:9",
               "https_proxy": "http://127.0.0.1:9", "all_proxy": "http://127.0.0.1:9", "NO_PROXY": ""}


def _resolve_ligand_bond_orders(src: Path, lig: dict, pdb_path: Path, rep) -> dict:
    """Try to establish bond orders for an extracted ligand (dictionary only; never guesses). Failure is not
    an error for receptor prep: the ligand is simply reported as needing a template."""
    from . import ligand_templates as lt
    from rdkit import Chem
    try:
        r = lt.ligand_from_structure(src, lig["resname"], lig["chain"], lig["seq"])
    except HelperError as e:
        rep.warn(f"ligand {lig['chain']}:{lig['resname']}{lig['seq']}: bond orders unknown ({e.message[:100]}...)")
        return {"status": "needs_template", "message": e.message, "needs": getattr(e, "needs", "smiles")}
    sdf = pdb_path.with_suffix(".sdf")
    w = Chem.SDWriter(str(sdf)); w.write(r.mol); w.close()
    return {"status": "resolved", "source": r.source, "sdf": str(sdf), "warnings": r.warnings}


def _run(cmd, cwd, timeout=None):
    """Run a helper tool. Stops it (whole process group) when the job is cancelled or exceeds the time limit."""
    import os
    import signal
    import time

    from . import cancel, limits
    from .errors import Cancelled
    cancel.check()
    limit = limits.STAGE_TIMEOUT_S if timeout is None else timeout
    proc = subprocess.Popen(cmd, cwd=cwd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True,
                            env={**os.environ, **OFFLINE_ENV, **resources.thread_env()}, start_new_session=True)
    deadline = time.time() + limit

    def stop():
        try:
            os.killpg(proc.pid, signal.SIGKILL)
        except (ProcessLookupError, PermissionError):
            proc.kill()
        proc.communicate()

    while True:
        try:
            out, err = proc.communicate(timeout=0.5)
            return subprocess.CompletedProcess(cmd, proc.returncode, out, err)
        except subprocess.TimeoutExpired:
            tok = cancel.current()
            if tok is not None and tok.cancelled:
                stop()
                raise Cancelled("job cancelled")
            if time.time() > deadline:
                stop()
                raise ToolFailure(f"{Path(str(cmd[0])).name} exceeded the {limit:.0f} s time limit and was stopped "
                                  "(very large structure? select fewer chains)")


def _pdbfixer(pdb_text, opts, rep):
    """Returns fixed PDB text."""
    from openmm.app import PDBFile
    from pdbfixer import PDBFixer
    fixer = PDBFixer(pdbfile=io.StringIO(pdb_text))
    fixer.findMissingResidues()
    chains = list(fixer.topology.chains())
    internal = {}
    for (ci, idx), seq in fixer.missingResidues.items():
        n = len(list(chains[ci].residues()))
        if 0 < idx < n:
            internal[(ci, idx)] = seq
    for (ci, idx), seq in internal.items():
        res = list(chains[ci].residues())[idx]
        rep.warn(f"missing internal loop of {len(seq)} residue(s) before chain "
                 f"{chains[ci].id} residue {res.name}{res.id}")
    if opts.model_loops_max > 0:
        fixer.missingResidues = {k: v for k, v in internal.items() if len(v) <= opts.model_loops_max}
        for (ci, idx), seq in fixer.missingResidues.items():
            rep.change("loops", f"{chains[ci].id}:{idx}", f"modelled {len(seq)}-residue loop",
                       "user requested short-loop modelling (speculative)", "PDBFixer")
    else:
        fixer.missingResidues = {}
    fixer.findNonstandardResidues()
    if fixer.nonstandardResidues:
        if opts.nonstandard == "map":
            for res, repl in fixer.nonstandardResidues:
                rep.change("nonstandard", f"{res.chain.id}:{res.name}{res.id}", f"mapped to {repl}",
                           "non-standard residue mapped to parent", "PDBFixer")
            fixer.replaceNonstandardResidues()
    if opts.rebuild_missing_atoms == "pdbfixer":
        fixer.findMissingAtoms()
        for res, atoms in fixer.missingAtoms.items():
            rep.change("missing_atoms", f"{res.chain.id}:{res.name}{res.id}",
                       "added " + ",".join(a.name for a in atoms), "incomplete side chain", "PDBFixer")
        # Terminal atoms (OXT) are left to PDB2PQR: PDBFixer's OXT placement can be badly off when the
        # carboxylate O is displaced (seen on real entries, e.g. 1CA2), which breaks Meeko's bonding.
        fixer.missingTerminals = {}
        fixer.addMissingAtoms(seed=0)
    out = io.StringIO()
    PDBFile.writeFile(fixer.topology, fixer.positions, out, keepIds=True)
    return out.getvalue()


def _pdb2pqr(pdb_path: Path, opts, rep, workdir: Path) -> Path:
    out = workdir / "protonated.pdb"
    cmd = [sys.executable, "-m", "pdb2pqr", "--ff", "AMBER", "--ffout", "AMBER", "--keep-chain",
           "--with-ph", str(opts.ph), "--titration-state-method", "propka",
           "--pdb-output", str(out)]
    if not opts.optimize_flips:
        cmd.append("--noopt")
    cmd += [str(pdb_path), str(workdir / "protonated.pqr")]
    r = _run(cmd, workdir)
    if r.returncode != 0 or not out.exists():
        m = re.search(r"atom count (\d+) exceeds (\d+)", r.stdout + r.stderr)
        if m:
            # Not a crash: the protonated model does not fit the PDB file format. No fallback can help.
            raise FormatLimit(f"the protonated model would have {int(m.group(1)):,} atoms, more than the PDB file format "
                              f"can hold ({int(m.group(2)):,}). Select fewer chains (give a reference ligand so only "
                              "the chains around the site are prepared, or use 'delete chain ...') and try again",
                              log=r.stdout + r.stderr)
        raise ToolFailure("PDB2PQR failed", log=r.stdout + r.stderr)
    log = r.stdout + r.stderr
    for m in re.finditer(r"Unable to debump (\w+) (\w) (\d+)", log):
        rep.warn(f"steric clash remained after hydrogen placement at {m.group(2)}:{m.group(1)}{m.group(3)}")
    rep.tools.append(record("pdb2pqr", _ver("pdb2pqr"), {"ph": opts.ph, "flips": opts.optimize_flips}))
    rep.tools.append(record("propka", _ver("propka"), {"ph": opts.ph}))
    return out


def _ver(pkg):
    from importlib.metadata import version
    try:
        return version(pkg)
    except Exception:  # noqa: BLE001
        return "unknown"


def _heavy_atoms(pdb_text):
    out = {}
    for line in pdb_text.splitlines():
        if line.startswith("ATOM") and line[76:78].strip() != "H" and not line[12:16].strip().startswith("H"):
            out.setdefault((line[21], line[22:27].strip()), set()).add(line[12:16].strip())
    return out


def _residue_names(pdb_text):
    out = {}
    for line in pdb_text.splitlines():
        if line.startswith("ATOM"):
            out[(line[21], line[22:27].strip())] = line[17:20].strip()
    return out


MAX_HEAVY_SHIFT = 0.5       # A: restrained minimisation must not move heavy atoms further than this
BOND_RANGE = {("N", "CA"): (1.35, 1.60), ("CA", "C"): (1.40, 1.65)}


def _topology_without_gap_bonds(pdb, app, unit):
    """OpenMM bonds consecutive residues blindly; across a chain break that pulls the gap shut."""
    top = app.Topology()
    chains = {}
    atoms = {}
    for chain in pdb.topology.chains():
        c = top.addChain(chain.id)
        for res in chain.residues():
            r = top.addResidue(res.name, c, res.id, res.insertionCode)
            for a in res.atoms():
                atoms[a.index] = top.addAtom(a.name, a.element, r, a.id)
    pos = pdb.positions
    for a, b in pdb.topology.bonds():
        if a.residue is not b.residue:
            d = sum((pos[a.index][k] - pos[b.index][k]).value_in_unit(unit.angstrom) ** 2 for k in range(3)) ** 0.5
            if d > 2.0:
                continue
        top.addBond(atoms[a.index], atoms[b.index])
    return top


_PARENT = {"HID": "HIS", "HIE": "HIS", "HIP": "HIS", "ASH": "ASP", "GLH": "GLU", "LYN": "LYS",
           "CYX": "CYS", "CYM": "CYS"}


def _parent_names(pdb_text: str) -> str:
    out = []
    for l in pdb_text.splitlines():
        if l.startswith(("ATOM", "HETATM")) and l[17:20] in _PARENT:
            l = l[:17] + _PARENT[l[17:20]] + l[20:]
        out.append(l)
    return "\n".join(out) + "\n"


def _geometry_problem(before: str, after: str) -> str:
    """Return a description if `after` is not a faithful, physically sane version of `before`."""
    def heavy(text):
        out = {}
        for l in text.splitlines():
            if l.startswith("ATOM") and not l[12:16].strip().startswith("H") and l[76:78].strip() != "H":
                out[(l[21], l[22:27].strip(), l[12:16].strip())] = (float(l[30:38]), float(l[38:46]), float(l[46:54]))
        return out
    b, a = heavy(before), heavy(after)
    shifts = {k: math.dist(b[k], a[k]) for k in b if k in a}
    worst_key = max(shifts, key=shifts.get) if shifts else None
    worst = shifts[worst_key] if worst_key else 0.0
    if worst > MAX_HEAVY_SHIFT:
        return f"heavy atom {worst_key[0]}:{worst_key[1]}.{worst_key[2]} moved {worst:.2f} A (limit {MAX_HEAVY_SHIFT})"
    for (ch, num), _ in {(k[0], k[1]): 0 for k in a}.items():
        for (n1, n2), (lo, hi) in BOND_RANGE.items():
            p, q = a.get((ch, num, n1)), a.get((ch, num, n2))
            if p and q and not lo <= math.dist(p, q) <= hi:
                return f"implausible {n1}-{n2} bond ({math.dist(p, q):.2f} A) at {ch}:{num}"
    return ""


PDB_FORMAT_ATOM_LIMIT = 99_999          # PDB files cannot number more atoms
HYDROGEN_FACTOR = 1.90                  # atoms after protonation per heavy atom, lower bound (6RD4: 2.02)
PDBFIXER_FALLBACK_MAX_RESIDUES = 1500   # PDBFixer needs minutes beyond this; fail with a clear message instead
MINIMIZE_ITERATIONS = 500
MINIMIZE_CHUNK = 25                 # iterations between checks
MINIMIZE_BUDGET_S = 240             # wall-clock cap for one minimisation (a 1,350-residue run took 490 s)


def _minimize_bounded(sim, pdb, topology, unit) -> int:
    """Minimise in small chunks; stop at once if a heavy atom drifts past the limit (the result would be rejected
    anyway; 9QGG spent 490 s to reach that verdict) or when the time budget is used up. Returns iterations done."""
    import time

    import numpy as np
    heavy = [a.index for a in topology.atoms() if a.element is not None and a.element.symbol != "H"]
    start = np.array(pdb.positions.value_in_unit(unit.nanometer))[heavy]
    t0, done = time.time(), 0
    while done < MINIMIZE_ITERATIONS:
        cancel.check()
        sim.minimizeEnergy(maxIterations=MINIMIZE_CHUNK)
        done += MINIMIZE_CHUNK
        now = np.array(sim.context.getState(getPositions=True).getPositions().value_in_unit(unit.nanometer))[heavy]
        shift = float(np.sqrt(((now - start) ** 2).sum(axis=1)).max()) * 10
        if shift > MAX_HEAVY_SHIFT:
            raise ValueError(f"result rejected, a heavy atom moved {shift:.2f} A within {done} iterations "
                             f"(limit {MAX_HEAVY_SHIFT}); stopped early")
        if time.time() - t0 > MINIMIZE_BUDGET_S:
            break
    return done


def _minimize(pdb_in: Path, pdb_out: Path, rep, k: float = RESTRAINT_K) -> bool:
    try:
        import openmm as mm
        from openmm import app, unit
        # OpenMM only bonds residue names it knows (not LYN, for one): hand it the parent names and let the
        # force-field templates tell variants apart by their hydrogens. Our names are restored on output.
        std = io.StringIO(_parent_names(pdb_in.read_text()))
        pdb = app.PDBFile(std)
        topology = _topology_without_gap_bonds(pdb, app, unit)
        ff = app.ForceField("amber14-all.xml", "implicit/gbn2.xml")
        # ignoreExternalBonds: residues at chain breaks have no neighbour and must still match.
        # Without external bonds OpenMM cannot choose between CYX (disulfide) and CYM (thiolate) for a Cys with
        # no HG, so pin by geometry: another SG within 2.5 A means a disulfide.
        sg = {r: pdb.positions[a.index] for r in topology.residues() for a in r.atoms() if a.name == "SG"}
        pinned = {}
        for r in topology.residues():
            if r.name == "CYS" and not any(a.name == "HG" for a in r.atoms()) and r in sg:
                near = any(o is not r and sum((sg[r][k] - sg[o][k]).value_in_unit(unit.angstrom) ** 2
                                              for k in range(3)) ** 0.5 < 2.5 for o in sg)
                pinned[r] = "CYX" if near else "CYM"
        system = ff.createSystem(topology, nonbondedMethod=app.CutoffNonPeriodic,
                                 nonbondedCutoff=1.5 * unit.nanometer, constraints=None,
                                 ignoreExternalBonds=True, residueTemplates=pinned)
        rest = mm.CustomExternalForce("0.5*k*((x-x0)^2+(y-y0)^2+(z-z0)^2)")
        rest.addGlobalParameter("k", k * unit.kilojoule_per_mole / unit.nanometer**2)
        for p in ("x0", "y0", "z0"):
            rest.addPerParticleParameter(p)
        for atom in pdb.topology.atoms():
            if atom.element.symbol != "H":
                rest.addParticle(atom.index, pdb.positions[atom.index].value_in_unit(unit.nanometer))
        system.addForce(rest)
        sim = app.Simulation(topology, system, mm.VerletIntegrator(0.001),
                             mm.Platform.getPlatformByName("CPU"), {"Threads": str(resources.cores())})
        sim.context.setPositions(pdb.positions)
        _minimize_bounded(sim, pdb, topology, unit)
        buf = io.StringIO()
        final_pos = sim.context.getState(getPositions=True).getPositions()
        app.PDBFile.writeFile(topology, final_pos, buf, keepIds=True)
        # OpenMM standardises HID/HIE/HIP, ASH, GLH, LYN, CYX... to their parent names; restore ours so
        # the protonation states survive into Meeko.
        names = _residue_names(pdb_in.read_text())
        lines = []
        for line in buf.getvalue().splitlines():
            if line.startswith(("ATOM", "HETATM")):
                nm = names.get((line[21], line[22:27].strip()))
                if nm:
                    line = line[:17] + nm.ljust(3) + line[20:]
            lines.append(line)
        text = "\n".join(lines) + "\n"
        problem = _geometry_problem(pdb_in.read_text(), text)
        if problem:
            raise ValueError(f"result rejected, {problem}")
        pdb_out.write_text(text)
        rep.change("minimization", "all", "restrained minimisation of hydrogens and side-chain strain",
                   f"heavy atoms restrained ({k / 4184:g} kcal/mol/A^2)", "OpenMM")
        rep.tools.append(record("openmm", _ver("openmm"), {"ff": "amber14+gbn2", "restraint_k_kj_mol_nm2": k}))
        return True
    except Exception as e:  # noqa: BLE001
        rep.warn(f"minimisation skipped: {e}")
        return False


def _meeko_templates() -> set:
    import json
    import meeko
    path = Path(meeko.__file__).parent / "data" / "residue_chem_templates.json"
    return set(json.loads(path.read_text())["residue_templates"])


def _report_pending_as_dropped(rep, dropped: dict, reason: str = "no bundled Meeko template") -> None:
    for (ch, num, name), n in dropped.items():
        rep.warn(f"{ch}:{name}{num} left out of the PDBQT: {reason} ({n} atoms)")
        rep.dropped.append({"residue": f"{ch}:{name}{num}", "reason": reason})


def _drop_untemplated(src: Path, dst: Path, rep, defer: bool = False) -> Path:
    """Remove residues Meeko has no bundled template for. Without this Meeko tries to download the
    definition from RCSB, which breaks offline use and the local-first rule (plan principle 5)."""
    templates = _meeko_templates()
    dropped, lines = {}, []
    text = src.read_text()
    has_hg = {(l[21], l[22:27].strip()) for l in text.splitlines()
              if l.startswith("ATOM") and l[12:16].strip() == "HG" and l[17:20] == "CYS"}
    renamed = set()
    for line in text.splitlines():
        if line.startswith(("ATOM", "HETATM")):
            name = line[17:20].strip()
            key = (line[21], line[22:27].strip())
            if name == "CYM" or (name == "CYS" and key not in has_hg):
                # Meeko has no thiolate template; CYX is the HG-less cysteine. The -1 charge is not represented.
                renamed.add(key)
                line = line[:17] + "CYX" + line[20:]
                name = "CYX"
            if name not in templates:
                dropped.setdefault((line[21], line[22:27].strip(), name), 0)
                dropped[(line[21], line[22:27].strip(), name)] += 1
                if defer:
                    rep.pending_het.setdefault((line[21], line[22:27].strip(), name), []).append(line)
                continue
        lines.append(line)
    for ch, num in sorted(renamed):
        rep.warn(f"{ch}:CYS{num} has no HG (thiolate or unbridged): represented as CYX in the PDBQT, "
                 "so its negative charge is not modelled")
    if not defer:
        _report_pending_as_dropped(rep, dropped)
    if not dropped and not renamed:
        return src
    dst.write_text("\n".join(lines) + "\n")
    return dst


MAX_VALENCE = {"H": 1, "C": 4, "N": 4, "O": 2, "S": 2}     # bonds a perceived-by-distance atom may have


def _malformed_residues(pdb_text: str) -> dict:
    """Residues whose own atoms sit so close that Meeko's distance rule gives an atom an impossible
    valence (e.g. a pentavalent carbon). That is a defect in the input coordinates, not something the
    pipeline can repair. Returns {(chain, resnum): (resname, "A.x - B.y (d A)")}."""
    from meeko.utils.covalent_radius_table import covalent_radius
    out = {}
    st = gemmi.read_pdb_string(pdb_text)
    for ch in st[0]:
        for res in ch:
            atoms = [a for a in res if a.element.name in covalent_radius]
            nbrs = {i: [] for i in range(len(atoms))}
            for i, a in enumerate(atoms):
                for j in range(i + 1, len(atoms)):
                    b = atoms[j]
                    d = a.pos.dist(b.pos)
                    if d < BOND_ALLOWANCE * (covalent_radius[a.element.name] + covalent_radius[b.element.name]):
                        nbrs[i].append((j, d))
                        nbrs[j].append((i, d))
            for i, a in enumerate(atoms):
                lim = MAX_VALENCE.get(a.element.name)
                if lim is not None and len(nbrs[i]) > lim:
                    j, d = min(nbrs[i], key=lambda t: t[1])
                    out[(ch.name, f"{res.seqid.num}{res.seqid.icode.strip()}")] = (
                        res.name, f"{a.name} has {len(nbrs[i])} bonds, closest {atoms[j].name} at {d:.2f} A")
                    break
    return out


def _drop_malformed(src: Path, dst: Path, rep) -> Path:
    """Leave out residues with impossible internal geometry, by name, instead of letting Meeko fail with
    an RDKit valence error. The result is reported as degraded (the receptor has a hole)."""
    text = src.read_text()
    bad = _malformed_residues(text)
    if not bad:
        return src
    for (ch, num), (name, why) in sorted(bad.items()):
        rep.warn(f"{ch}:{name}{num} left out of the PDBQT: impossible geometry in the input structure ({why})")
        rep.dropped.append({"residue": f"{ch}:{name}{num}", "reason": f"impossible geometry in input ({why})"})
    keep = [l for l in text.splitlines()
            if not (l.startswith(("ATOM", "HETATM")) and (l[21], l[22:27].strip()) in bad)]
    dst.write_text("\n".join(keep) + "\n")
    return dst


CROSSCHECK_MAX_RESIDUES = 2500      # the optional Meeko cross-check is skipped above this size


def _write_pdbqt(pdb_path: Path, base: Path, opts, rep):
    """Table writer first (it cannot be broken by odd geometry), Meeko as fallback or cross-check."""
    from . import recwriter
    mode = opts.pdbqt_writer
    if mode not in ("auto", "table", "meeko"):
        raise InvalidInput(f"pdbqt_writer must be auto, table or meeko, not {mode!r}")
    if opts.flex_residues:
        if mode == "table":
            raise InvalidInput("flexible residues need Meeko (pdbqt_writer='table' writes a rigid receptor only)")
        mode = "meeko"                       # the table writer has no torsion trees
    if mode == "meeko":
        rep.pdbqt_info = {"writer": "meeko", "crosscheck": "not applicable"}
        _settle_pending_het(rep, None)
        return _meeko(pdb_path, base, opts, rep)

    out = Path(f"{base}.pdbqt")
    r = recwriter.write_pdbqt(pdb_path.read_text(), out)
    if r["dropped"] and mode == "auto":
        # The table lacks rules for something (modified residue, nucleic acid, ...). Meeko may know it.
        trial = _Report()
        alt_base = base.parent / "meeko_full" / "receptor"
        alt_base.parent.mkdir(exist_ok=True)
        try:
            alt, _ = _meeko(pdb_path, alt_base, opts, trial)
            if len(trial.dropped) < len(r["dropped"]):
                import shutil
                shutil.copy(alt, out)
                rep.changes += trial.changes
                rep.warnings += [w for w in trial.warnings if w not in rep.warnings]
                rep.dropped += trial.dropped
                rep.tools += trial.tools
                rep.warn("the built-in table writer has no rules for: "
                         + ", ".join(d["residue"] for d in r["dropped"][:5]) + "; Meeko wrote the PDBQT instead")
                rep.pdbqt_info = {"writer": "meeko", "crosscheck": "not applicable (table writer incomplete)"}
                _settle_pending_het(rep, None)
                return out, None
        except ToolFailure:
            pass
    for d in r["dropped"]:
        rep.warn(f"residue {d['residue']} left out of the PDBQT: {d['reason']}")
        rep.dropped.append(d)
    rep.tools.append(record("dockchem.recwriter", __version__, {"atoms": r["atoms"], "merged_hydrogens": r["merged_hydrogens"],
                                                                 "charges": "not written (Vina ignores them)"}))
    rep.pdbqt_info = {"writer": "table", "atoms": r["atoms"]}
    _settle_pending_het(rep, out, r["atoms"])
    if mode == "auto" and opts.crosscheck:
        n_res = len({(l[21], l[22:27]) for l in pdb_path.read_text().splitlines() if l.startswith(("ATOM", "HETATM"))})
        if n_res > CROSSCHECK_MAX_RESIDUES:
            rep.pdbqt_info["crosscheck"] = {"result": "skipped", "reason": f"{n_res} residues (limit {CROSSCHECK_MAX_RESIDUES})"}
        else:
            rep.pdbqt_info["crosscheck"] = _crosscheck(pdb_path, out, base, opts, rep)
    return out, None


def _settle_pending_het(rep, out: Path | None, serial: int = 0) -> None:
    """Kept hetero groups Meeko has no template for: type them through the ligand route and append them to the
    table-written PDBQT; anything that cannot be typed is reported as left out (status degraded)."""
    from . import recwriter
    pending, rep.pending_het = rep.pending_het, {}
    unsupported = {}
    typed = []
    for (ch, num, name), lines in pending.items():
        if out is None:
            unsupported[(ch, num, name)] = (len(lines), "no bundled Meeko template")
            continue
        try:
            seq = int(re.match(r"-?\d+", num).group())
            if len(lines) == 1 and lines[0][76:78].strip().capitalize() in AD4_METALS:    # single ion with an AD4 type
                el = lines[0][76:78].strip().capitalize()
                serial = recwriter.append_hetero(out, ch, name, seq, [(lines[0][12:16].strip() or el, (
                    float(lines[0][30:38]), float(lines[0][38:46]), float(lines[0][46:54])), AD4_METALS[el])], serial)
                typed.append(f"{ch}:{name}{num}")
                continue
            atoms, warns = recwriter.type_hetero_group("\n".join(lines) + "\nEND\n", name, ch, seq)
            serial = recwriter.append_hetero(out, ch, name, seq, atoms, serial)
            typed.append(f"{ch}:{name}{num}")
            for w in warns:
                rep.warn(f"{ch}:{name}{num}: {w}")
        except HelperError as e:
            unsupported[(ch, num, name)] = (len(lines), e.message.split(". ")[0][:150])
    for (ch, num, name), (n, why) in unsupported.items():
        rep.warn(f"{ch}:{name}{num} left out of the PDBQT: {why} ({n} atoms)")
        rep.dropped.append({"residue": f"{ch}:{name}{num}", "reason": why})
    if typed:
        rep.change("hetero_groups", ", ".join(typed), "added to the receptor PDBQT",
                   "kept groups typed from the bundled wwPDB dictionary", "dockchem.recwriter")


def _crosscheck(pdb_path: Path, table_pdbqt: Path, base: Path, opts, rep) -> dict:
    """Run Meeko independently and compare AutoDock types per residue. Never fatal."""
    import shutil
    work = base.parent / "meeko_crosscheck"
    work.mkdir(exist_ok=True)
    try:
        throwaway = _Report()
        meeko_out, _ = _meeko(pdb_path, work / "receptor", opts, throwaway)
        result = compare_pdbqt(table_pdbqt, meeko_out)
    except ToolFailure as e:
        result = {"result": "unavailable", "reason": e.message}
    except Exception as e:  # noqa: BLE001
        result = {"result": "unavailable", "reason": f"{type(e).__name__}: {e}"}
    finally:
        shutil.rmtree(work, ignore_errors=True)
    rep.tools.append(record("meeko", _ver("meeko"), {"role": "cross-check", "result": result["result"]}))
    if result["result"] == "disagree":
        rep.warn(f"Meeko and the built-in writer disagree on {result['n_residues_different']} residue(s), first: "
                 f"{result['examples'][0]}")
    return result


def compare_pdbqt(table_path: Path, meeko_path: Path) -> dict:
    """Per-residue comparison of AutoDock atom types (atom names differ: Meeko renames N-terminal H and water
    atoms). Residues Meeko left out are not differences; the report lists them separately."""
    import collections

    def per_residue(path):
        d = collections.defaultdict(collections.Counter)
        for l in Path(path).read_text().splitlines():
            if l.startswith(("ATOM", "HETATM")):
                d[(l[21], l[22:27].strip(), l[17:20].strip())][l[77:79].strip()] += 1
        return d
    a, b = per_residue(table_path), per_residue(meeko_path)
    shared = set(a) & set(b)
    conflicts, extra, hydrogen_only = [], [], []
    for k in sorted(shared):
        missing_in_table = b[k] - a[k]          # atom types Meeko wrote that the table did not
        only_table = a[k] - b[k]                # atoms the table wrote that Meeko left out (terminal atoms)
        label = f"{k[0]}:{k[2]}{k[1]}"
        if missing_in_table and set(missing_in_table) <= {"HD"} and not set(only_table) - {"HD"}:
            hydrogen_only.append(f"{label} polar H table {a[k]['HD']} vs Meeko {b[k]['HD']}")
        elif missing_in_table:
            conflicts.append(f"{label} table {dict(a[k])} vs Meeko {dict(b[k])}")
        elif only_table:
            extra.append(f"{label} +{dict(only_table)}")
    return {"result": "disagree" if conflicts else "agree", "n_residues_compared": len(shared),
            "n_residues_different": len(conflicts), "examples": conflicts[:10],
            "n_residues_hydrogen_count_differs": len(hydrogen_only), "hydrogen_examples": hydrogen_only[:10],
            "n_residues_table_has_extra_atoms": len(extra), "extra_atoms_examples": extra[:5],
            "only_in_table": sorted(f"{k[0]}:{k[2]}{k[1]}" for k in set(a) - set(b))[:20],
            "only_in_meeko": sorted(f"{k[0]}:{k[2]}{k[1]}" for k in set(b) - set(a))[:20]}


def _isolate_valence_failures(pdb: Path, attempt) -> dict:
    """{"A:46": "ASN"} for every residue that makes Meeko fail on its own with a valence error."""
    atoms = [l for l in pdb.read_text().splitlines() if l.startswith(("ATOM", "HETATM"))]
    order, names = [], {}
    for l in atoms:
        k = f"{l[21].strip()}:{l[22:27].strip()}"
        if k not in names:
            order.append(k)
            names[k] = l[17:20].strip()
    probe = pdb.with_name(pdb.stem + "_probe.pdb")

    def fails(keys):
        wanted = set(keys)
        probe.write_text("\n".join(l for l in atoms if f"{l[21].strip()}:{l[22:27].strip()}" in wanted) + "\nEND\n")
        r = attempt(probe)
        return r.returncode != 0 and "valence" in r.stdout + r.stderr

    def find(keys):
        if not fails(keys):
            return []
        if len(keys) == 1:
            return list(keys)
        h = len(keys) // 2
        return find(keys[:h]) + find(keys[h:])

    return {k: names[k] for k in find(order)}


def _meeko(pdb_path: Path, base: Path, opts, rep):
    """Run mk_prepare_receptor. Ladder: as-is -> without OXT on rejected C-termini -> delete (named)."""
    exe = Path(sys.executable).parent / "mk_prepare_receptor.py"
    exe = str(exe) if exe.exists() else "mk_prepare_receptor.py"
    flex = [a for f in opts.flex_residues for a in ("-f", f)]

    def attempt(path: Path, extra=()):
        return _run([exe, "-i", str(path), "-o", str(base), "-p", *flex, *extra], base.parent)

    def failing(res):
        return sorted(set(re.findall(r"No template matched for residue_key='([^']+)'", res.stdout + res.stderr)))

    used, deleted = pdb_path, False
    r = attempt(used)
    if r.returncode != 0 and "valence" in r.stdout + r.stderr:
        # RDKit rejected a residue's own geometry (e.g. a pentavalent carbon from distorted coordinates).
        # Find the offending residues by bisection, leave them out by name, and carry on.
        bad = _isolate_valence_failures(used, attempt)
        if bad:
            keep = [l for l in used.read_text().splitlines()
                    if not (l.startswith(("ATOM", "HETATM")) and f"{l[21].strip()}:{l[22:27].strip()}" in bad)]
            alt = used.with_name(used.stem + "_sane.pdb")
            alt.write_text("\n".join(keep) + "\n")
            for key, name in sorted(bad.items()):
                label = key.replace(":", f":{name}", 1)
                rep.warn(f"residue {label} left out of the PDBQT: its coordinates give an impossible "
                         "valence (distorted geometry in the input structure)")
                rep.dropped.append({"residue": label, "reason": "impossible valence (distorted input geometry)"})
            used = alt
            r = attempt(used)
    if r.returncode != 0:
        keys = set(failing(r))
        lines, stripped = [], set()
        for l in used.read_text().splitlines():
            key = f"{l[21].strip()}:{l[22:27].strip()}" if l.startswith(("ATOM", "HETATM")) else None
            if key in keys and l[12:16].strip() == "OXT":
                stripped.add(key)
                continue
            lines.append(l)
        if stripped:
            alt = pdb_path.with_name(pdb_path.stem + "_noOXT.pdb")
            alt.write_text("\n".join(lines) + "\n")
            r_alt = attempt(alt)
            if r_alt.returncode == 0:
                for key in sorted(stripped):
                    rep.warn(f"residue {key}: OXT omitted from the PDBQT because Meeko could not match "
                             "the C-terminal template with it")
                r, used = r_alt, alt
            else:
                r, used = r_alt, alt                      # keep the smaller failing set for the next step
    if r.returncode != 0:
        names = failing(r)
        r2 = attempt(used, ["-x"])        # last resort: delete residues that match no template
        if r2.returncode != 0:
            hint = " (hydrogen clash left by PDB2PQR?)" if "valence" in (r.stdout + r.stderr) else ""
            raise ToolFailure("Meeko receptor preparation failed" + hint, log=r.stdout + r.stderr)
        deleted = True
        for key in names:
            rep.warn(f"residue {key} left out of the PDBQT: no Meeko template matched (clashes or unusual atoms)")
            rep.dropped.append({"residue": key, "reason": "no Meeko template matched (clashes or unusual atoms)"})
        r = r2
    rep.tools.append(record("meeko", _ver("meeko"), {"flex": opts.flex_residues, "delete_bad_res": deleted}))
    if opts.flex_residues:
        rigid, flexf = Path(f"{base}_rigid.pdbqt"), Path(f"{base}_flex.pdbqt")
        if not rigid.exists():
            raise ToolFailure("Meeko did not produce rigid/flex PDBQT", log=r.stdout + r.stderr)
        return rigid, flexf
    out = Path(f"{base}.pdbqt")
    if not out.exists():
        raise ToolFailure("Meeko produced no PDBQT", log=r.stdout + r.stderr)
    return out, None


BOND_ALLOWANCE = 1.2       # Meeko perceives a bond when d < allowance * (r1 + r2) (meeko.polymer)


def _contact_pairs(pdb_text: str) -> list[tuple]:
    """Inter-residue atom pairs that Meeko would read as covalent bonds, using Meeko's own rule, as
    ((chain, resname, seq, icode, atom), (...), dist)."""
    from meeko.utils.covalent_radius_table import covalent_radius
    st = gemmi.read_pdb_string(pdb_text)
    ns = gemmi.NeighborSearch(st[0], st.cell, 5).populate()
    cs = gemmi.ContactSearch(2 * BOND_ALLOWANCE * max(covalent_radius.values()))
    cs.ignore = gemmi.ContactSearch.Ignore.SameResidue      # adjacent residues clash too (3SN6 THR98 HG1)
    out = []
    for r in cs.find_contacts(ns):
        a, b = r.partner1, r.partner2
        if a.atom.name == "SG" and b.atom.name == "SG":      # disulfide
            continue
        if {a.atom.name, b.atom.name} == {"C", "N"} and r.dist > 1.0:     # the peptide bond itself
            continue
        ea, eb = a.atom.element.name, b.atom.element.name       # table is keyed by element symbol
        if ea not in covalent_radius or eb not in covalent_radius:
            continue
        if r.dist >= BOND_ALLOWANCE * (covalent_radius[ea] + covalent_radius[eb]):
            continue
        ident = lambda p: (p.chain.name, p.residue.name, p.residue.seqid.num,
                           p.residue.seqid.icode.strip(), p.atom.name)
        out.append((ident(a), ident(b), r.dist))
    return out


def _close_contacts(pdb_text: str) -> list[str]:
    return [f"{a[0]}:{a[1]}{a[2]}.{a[4]} - {b[0]}:{b[1]}{b[2]}.{b[4]} ({d:.2f} A)"
            for a, b, d in _contact_pairs(pdb_text)]


def _strip_clashing_hydrogens(pdb_text: str, rep) -> tuple[str, int]:
    """Last resort: drop only the hydrogens that sit on top of other atoms, so Meeko can still match the
    residue (a deleted residue is far worse than a missing hydrogen). Every removal is reported."""
    drop = set()
    for a, b, _ in _contact_pairs(pdb_text):
        for p in (a, b):
            if p[4].startswith("H"):
                drop.add((p[0], f"{p[2]}{p[3]}", p[4]))
                break                       # one hydrogen per contact is enough
    if not drop:
        return pdb_text, 0
    kept = [l for l in pdb_text.splitlines()
            if not (l.startswith("ATOM") and (l[21], l[22:27].strip(), l[12:16].strip()) in drop)]
    for ch, num, name in sorted(drop):
        rep.change("hydrogen_repair", f"{ch}:{num}", f"removed {name}", "could not be repositioned without moving heavy atoms", "dockchem")
    return "\n".join(kept) + "\n", len(drop)


def _clashes(pdb_text, rep):
    st = gemmi.read_pdb_string(pdb_text)
    st.remove_hydrogens()
    ns = gemmi.NeighborSearch(st[0], st.cell, 5).populate()
    seen = set()
    cs = gemmi.ContactSearch(CLASH_CUTOFF)
    cs.ignore = gemmi.ContactSearch.Ignore.AdjacentResidues
    for r in cs.find_contacts(ns):
        if r.partner1.atom.name == "SG" and r.partner2.atom.name == "SG":      # disulfide bond, not a clash
            continue
        key = tuple(sorted((r.partner1.residue.name + str(r.partner1.residue.seqid.num),
                            r.partner2.residue.name + str(r.partner2.residue.seqid.num))))
        if key not in seen:
            seen.add(key)
            rep.warn(f"steric clash: {key[0]} / {key[1]} ({r.dist:.2f} A)")


def prepare_receptor(src: str | Path, outdir: str | Path,
                     opts: ReceptorPrepOptions | None = None) -> ReceptorResult:
    opts = opts or ReceptorPrepOptions()
    src, outdir = Path(src), Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    rep = _Report()
    from .limits import check_structure_text
    check_structure_text(src)
    try:
        st = gemmi.read_structure(str(src))
    except Exception as e:  # noqa: BLE001
        raise InvalidInput(f"cannot read structure: {e}") from e
    if len(st) == 0:
        raise InvalidInput("structure has no models")
    st.setup_entities()
    if not 1 <= opts.model <= len(st):
        raise InvalidInput(f"model {opts.model} not in structure (has {len(st)})")
    if len(st) > 1:
        n_models = len(st)
        for i in range(n_models - 1, -1, -1):
            if i != opts.model - 1:
                del st[i]
        rep.change("model", f"model {opts.model}", f"kept 1 of {n_models} models", "single model docked", "gemmi")
    model = st[0]
    assessment = assess_structure(st, 1, src.suffix.lower())
    n_h = sum(1 for ch in model for r in ch for a in r if a.element.is_hydrogen)
    atoms_before = sum(len(r) for ch in model for r in ch)
    if n_h:
        st.remove_hydrogens()
        rep.change("hydrogens", "all", f"removed {n_h} input hydrogens", "re-added consistently by PDB2PQR", "gemmi")
    model = st[0]

    # site residue occupancy / B-factor diagnostics need pre-collapse info
    partial_occ = {_rid(ch.name, r) for ch in model for r in ch if any(a.occ < 1.0 for a in r)}
    _collapse_altlocs(model, rep)

    if opts.edits:
        from .selection import EditSession
        sess = EditSession(st)
        for cmd in opts.edits:
            e = sess.apply(cmd)
            rep.change("selection", f"{e.n_residues} residue(s), {e.n_atoms} atom(s)", cmd,
                       "user edit", "dockchem.selection")
            for w in e.warnings:
                rep.warn(f"edit '{cmd}': {w}")
        st = sess.materialize()
        model = st[0]
        if not any(_is_polymer(r) for ch in model for r in ch):
            raise InvalidInput("the edits removed every protein residue")
        assessment = assess_structure(st, 1, ".cif")        # what is left is what gets docked into
    if assessment["verdict"] == "unusable":
        if not opts.allow_unusable:
            raise UnusableStructure("this structure cannot be docked into: " + assessment["summary"] +
                                    " (set allow_unusable to prepare it anyway)", assessment)
        rep.warn("structure assessed as unusable but prepared on request: " + assessment["summary"])
    opts, chain_map = _alias_chains(st, opts, rep)
    model = st[0]
    refs_all = _find_reference(model, opts.reference_ligand)
    refs = refs_all[:1]
    if len(refs_all) > 1:
        rep.warn(f"{len(refs_all)} copies of {opts.reference_ligand} found; site defined by "
                 f"{_rid(refs[0][0], refs[0][1])}; others: "
                 + ", ".join(_rid(c, r) for c, r in refs_all[1:]) + " (plan §8.2: pick one)")
    ref_pos = [a.pos for _, r in refs for a in r]
    if opts.reference_ligand and not refs:
        raise InvalidInput(f"reference ligand {opts.reference_ligand} not found"
                           + (" (was it removed by one of the edits?)" if opts.edits else ""))

    # chain selection
    ns = gemmi.NeighborSearch(model, st.cell, 6).populate()
    poly_chains = [ch.name for ch in model if any(_is_polymer(r) for r in ch)]
    if opts.chains is not None:
        chains = list(opts.chains)
        bad = [c for c in chains if c not in poly_chains]
        if bad:
            raise InvalidInput(f"chain(s) {bad} not found or not protein")
    elif ref_pos:
        touching = set()
        for p in ref_pos:
            for m in ns.find_atoms(p, "\0", radius=5.0):
                cra = m.to_cra(model)
                if _is_polymer(cra.residue):
                    touching.add(cra.chain.name)
        chains = [c for c in poly_chains if c in touching] or poly_chains
    else:
        chains = poly_chains
    dropped = [c for c in poly_chains if c not in chains]
    for c in dropped:
        rep.change("chains", c, "removed chain", "not within 5 A of reference ligand / not selected", "gemmi")
    if not chains:
        raise InvalidInput("no protein chains found")

    # classify heteroatom groups
    site_pos = ref_pos
    ref_keys = {(c, r.name, r.seqid.num, r.seqid.icode) for c, r in refs}
    ligands = []
    het_out: list[tuple[str, gemmi.Residue]] = []
    for ch in model:
        for res in ch:
            if _is_polymer(res):
                continue
            cls = hetatm.classify(res.name, _heavy(res), opts.hetatm_overrides)
            rid = _rid(ch.name, res)
            in_chain = ch.name in chains
            near = (not site_pos) or _min_dist(res, site_pos) <= opts.site_radius
            if cls == "water":
                if opts.water_keep_radius is not None and site_pos and in_chain and \
                        _min_dist(res, site_pos) <= opts.water_keep_radius:
                    het_out.append((ch.name, res))
                    rep.change("waters", rid, "kept", f"within {opts.water_keep_radius:g} A of ligand", "gemmi")
                continue
            is_ref = (ch.name, res.name, res.seqid.num, res.seqid.icode) in ref_keys
            if is_ref or cls == "ligand":
                # The reference ligand is always extracted, even if its name is on a cofactor list.
                if is_ref or not refs:
                    ligands.append({"chain": ch.name, "resname": res.name, "seq": res.seqid.num,
                                    "heavy_atoms": _heavy(res)})
                rep.change("hetatm", rid, "removed from receptor", "ligand of interest (kept as separate file)", "gemmi+hetatm")
            elif cls in ("metal", "cofactor", "glycan") and in_chain and (cls == "glycan" or near):
                het_out.append((ch.name, res))
                rep.change("hetatm", rid, f"kept ({cls})", "in site / part of structure", "gemmi+hetatm")
                if cls == "metal":
                    rep.warn(f"metal {rid} retained: Vina models metal interactions poorly"
                             + ("; consider AD4Zn for zinc" if res.name == "ZN" else ""))
                if not site_pos and cls != "glycan":
                    rep.warn(f"{rid} kept without a reference site (no ligand to define 'in site')")
            else:
                rep.change("hetatm", rid, f"removed ({cls})", "additive / outside site / unselected chain", "gemmi+hetatm")

    blockers = [(c, r, why) for c, r in het_out if (why := _unsupported_reason(r))
                and site_pos and _min_dist(r, site_pos) <= SITE_IMPACT_RADIUS]
    if blockers:
        names = ", ".join(f"{_rid(c, r)} ({why})" for c, r, why in blockers[:4])
        msg = (f"the binding site contains a group that this version cannot represent: {names}. Docking into a "
               "pocket missing it would give misleading results")
        if not opts.allow_unusable:
            raise UnusableStructure(msg + " (set allow_unusable to prepare it anyway)", {
                "verdict": "unusable", "summary": msg,
                "findings": [{"severity": "block", "code": "unsupported_site_cofactor", "message": msg,
                              "residues": [_rid(c, r) for c, r, _ in blockers]}], "stats": {}})
        rep.warn("prepared on request although " + msg)

    # write ligands before editing
    for lig in ligands:
        lig_st = gemmi.Structure()
        lm = gemmi.Model("1")
        lc = gemmi.Chain(lig["chain"])
        for ch in model:
            for res in ch:
                if ch.name == lig["chain"] and res.name == lig["resname"] and res.seqid.num == lig["seq"]:
                    lc.add_residue(res)
        lm.add_chain(lc)
        lig_st.add_model(lm)
        p = outdir / f"ligand_{lig['chain']}_{lig['resname']}{lig['seq']}.pdb"
        lig_st.write_pdb(str(p))
        lig["path"] = str(p)
        lig["bond_orders"] = _resolve_ligand_bond_orders(src, lig, p, rep)

    # build polymer-only structure for PDBFixer / PDB2PQR
    poly = gemmi.Structure()
    poly.spacegroup_hm, poly.cell, poly.name = st.spacegroup_hm, st.cell, st.name
    pm = gemmi.Model("1")
    for ch in model:
        if ch.name not in chains:
            continue
        pc = gemmi.Chain(ch.name)
        no_ring = {(r.seqid.num, r.seqid.icode) for r in truncated_segment_start_prolines(ch)}
        for res in ch:
            if _is_polymer(res):
                if res.name in hetatm.MODIFIED and opts.nonstandard == "remove":
                    rep.change("nonstandard", _rid(ch.name, res), "removed", "user chose to remove", "gemmi")
                    continue
                stub = _is_backbone_stub(res) if gemmi.find_tabulated_residue(res.name) else None
                if not stub and (res.seqid.num, res.seqid.icode) in no_ring:
                    stub = "proline at a chain start without its ring atoms CG/CD (PDB2PQR crashes on it)"
                if stub:
                    rep.warn(f"{_rid(ch.name, res)} left out: {stub}")
                    rep.dropped.append({"residue": _rid(ch.name, res), "reason": stub})
                    continue
                pc.add_residue(res)
                pc[len(pc) - 1].het_flag = "A"
        pm.add_chain(pc)
    poly.add_model(pm)
    poly.setup_entities()

    site_res_ids = set()
    if site_pos:
        for ch in pm:
            for r in ch:
                if _min_dist(r, site_pos) <= opts.site_radius:
                    site_res_ids.add(_rid(ch.name, r))
    if opts.alphafold:
        _alphafold(pm, site_res_ids, opts, rep)
    for rid in sorted(partial_occ & site_res_ids):
        rep.warn(f"occupancy < 1 in site: {rid}")
    if site_res_ids:
        bs = [a.b_iso for ch in pm for r in ch for a in r]
        if len(bs) > 5:
            mean = sum(bs) / len(bs)
            sd = (sum((b - mean) ** 2 for b in bs) / len(bs)) ** 0.5
            for ch in pm:
                for r in ch:
                    if _rid(ch.name, r) in site_res_ids and sd > 0 and any(a.b_iso > mean + 3 * sd for a in r):
                        rep.warn(f"B-factor outlier in site: {_rid(ch.name, r)}")
    for ch_name, a, b in _chain_breaks(pm):
        near = site_pos and (_min_dist(a, site_pos) <= LOOP_NEAR_SITE or _min_dist(b, site_pos) <= LOOP_NEAR_SITE)
        rep.warn(f"chain break between {_rid(ch_name, a)} and {_rid(ch_name, b)}"
                 + (" within 10 A of the site: the pocket may be incomplete" if near else ""))

    with tempfile.TemporaryDirectory() as td:
        td = Path(td)
        cancel.check()
        fixed = _pdbfixer(poly.make_pdb_string(), opts, rep)
        stagecheck.check_moved_atoms("pdbfixer", poly.make_pdb_string(), fixed, limit=stagecheck.PDBFIXER_MOVE_LIMIT)
        (td / "fixed.pdb").write_text(fixed)
        rep.tools.append(record("pdbfixer", _ver("pdbfixer"), {"rebuild_missing_atoms": opts.rebuild_missing_atoms}))
        before = _residue_names(fixed)
        n_heavy = sum(1 for l in fixed.splitlines() if l.startswith(("ATOM", "HETATM")) and l[76:78].strip() not in ("H", "D"))
        if n_heavy * HYDROGEN_FACTOR > PDB_FORMAT_ATOM_LIMIT:
            raise FormatLimit(f"{n_heavy:,} heavy atoms would give a protonated model of about {int(n_heavy * 2):,} atoms, more "
                              f"than the PDB file format can hold ({PDB_FORMAT_ATOM_LIMIT:,}). Select fewer chains (give a "
                              "reference ligand so only the chains around the site are prepared, or use 'delete chain ...') "
                              "and try again")
        bad_dec = {k: v for k, v in opts.flip_decisions.items() if str(v).lower() not in ("accept", "reject")}
        if bad_dec:
            raise InvalidInput(f"flip_decisions must be 'accept' or 'reject': {bad_dec}")
        fixed, applied = protonation.apply_overrides(fixed, opts.protonation_overrides)
        (td / "fixed.pdb").write_text(fixed)
        user_set = {tuple(a["key"].split(":", 1)) for a in applied}
        try:
            prot = _pdb2pqr(td / "fixed.pdb", opts, rep, td)
        except ToolFailure as first:
            if opts.rebuild_missing_atoms != "pdb2pqr":
                raise
            n_res = len({(l[21], l[22:27]) for l in fixed.splitlines() if l.startswith("ATOM")})
            if n_res > PDBFIXER_FALLBACK_MAX_RESIDUES:
                raise ToolFailure(f"PDB2PQR failed and the structure ({n_res} residues) is too large for the slower "
                                  f"PDBFixer fallback (limit {PDBFIXER_FALLBACK_MAX_RESIDUES}); select the chains you need "
                                  "(a reference ligand or 'delete chain ...') and try again", log=first.log) from first
            # PDB2PQR dies on some truncated residues (e.g. an N-terminal Pro without its ring). Let PDBFixer
            # complete the heavy atoms and give PDB2PQR whole residues.
            import dataclasses
            alt = dataclasses.replace(opts, rebuild_missing_atoms="pdbfixer")
            fixed = _pdbfixer(poly.make_pdb_string(), alt, rep)
            fixed, applied = protonation.apply_overrides(fixed, opts.protonation_overrides)
            (td / "fixed.pdb").write_text(fixed)
            try:
                prot = _pdb2pqr(td / "fixed.pdb", alt, rep, td)
            except ToolFailure:
                raise first
            rep.warn("PDB2PQR could not process the input as given; missing heavy atoms were added by PDBFixer "
                     "instead and PDB2PQR then protonated the completed residues")
        prot_text = prot.read_text()
        protonation.verify_overrides(prot_text, applied)
        flips = protonation.detect_flips(fixed, prot_text) if opts.optimize_flips else []
        flipped = {f["key"] for f in flips}
        doomed = {tuple(k.split(":", 1)) for k in assessment["stats"].get("distorted_keys", [])}   # dropped later anyway
        stagecheck.check_moved_atoms("pdb2pqr", fixed, prot_text,
                                     {tuple(f["key"].split(":", 1)) for f in flips} | doomed)
        bare = stagecheck.residues_without_hydrogens(prot_text)
        if bare:
            # Residues PDB2PQR has no template for (MLY, SEC, ASX, UNK, ...) come back without hydrogens. Leave them
            # out by name, like other unsupported residues, rather than failing the whole entry.
            prot_text = stagecheck.drop_residues(prot_text, bare)
            prot.write_text(prot_text)
            for c, n, res in bare:
                rep.warn(f"{c}:{res}{n} left out: PDB2PQR has no protonation template for {res}")
                rep.dropped.append({"residue": f"{c}:{res}{n}", "reason": f"no protonation template for {res}"})
        stagecheck.check_hydrogens("pdb2pqr", prot_text)
        reject = []
        for k, v in opts.flip_decisions.items():
            c, n = protonation.parse_key(k)
            if str(v).lower() == "reject":
                if f"{c}:{n}" in flipped:
                    reject.append(f"{c}:{n}")
                else:
                    rep.warn(f"flip decision {k}: PDB2PQR did not flip this residue, nothing to reject")
        if reject:
            prot_text = protonation.undo_flips(prot_text, reject)
            prot.write_text(prot_text)
        for f in flips:
            c, n = f["key"].split(":", 1)
            kept = f["key"] in reject
            rep.change("flip", f"{c}:{before.get(('%s' % c, n), f['resname'])}{n}",
                       "flip rejected, crystal orientation restored" if kept else "side chain flipped",
                       "user decision" if kept else "hydrogen-bond optimisation (PDB2PQR); "
                       "set flip_decisions to reject", "user" if kept else "PDB2PQR")
            if kept:
                rep.warn(f"flip at {f['key']} was rejected: the crystal orientation may form worse "
                         "hydrogen bonds")
        if not opts.optimize_flips and opts.flip_decisions:
            rep.warn("flip_decisions ignored because optimize_flips is off (no flips are made)")
        if opts.rebuild_missing_atoms == "pdb2pqr":
            old_h, new_h = _heavy_atoms(fixed), _heavy_atoms(prot_text)
            for key, names in new_h.items():
                added = sorted(names - old_h.get(key, set()))
                if added:
                    rep.change("missing_atoms", f"{key[0]}:{before.get(key, '?')}{key[1]}",
                               "added " + ",".join(added), "missing heavy/terminal atoms", "PDB2PQR")
        for key, new in _residue_names(prot_text).items():
            old = before.get(key)
            if old and old != new:
                mine = key in user_set
                rep.change("protonation", f"{key[0]}:{old}{key[1]}", f"{old} -> {new}",
                           "user override" if mine else
                           f"PROPKA pKa at pH {opts.ph:g}" if new != "HID" and new != "HIE" else
                           "histidine tautomer from H-bond network",
                           "user" if mine else "PDB2PQR+PROPKA")
                if new in UNUSUAL_STATES:
                    rep.warn(f"unusual protonation {key[0]}:{new}{key[1]} ({UNUSUAL_STATES[new]})"
                             + (" - set by user override" if mine else ""))
        cancel.check()
        final = outdir / "receptor_prepared.pdb"
        if opts.minimize:
            if not _minimize(prot, final, rep):
                final.write_text(prot_text)
        else:
            final.write_text(prot_text)
            bad = _close_contacts(prot_text)
            if bad:
                # PDB2PQR can leave hydrogens on top of acceptors (salt bridges); Meeko then reads them as
                # bonds. Repair with restrained minimisation (heavy atoms held) rather than failing or
                # silently deleting residues. Deviates from "minimisation off by default" (plan §5.1).
                rep.warn(f"{len(bad)} implausibly close contact(s) after hydrogen placement, first: {bad[0]}")
                if _minimize(prot, final, rep, REPAIR_K):
                    rep.changes[-1].update(stage="hydrogen_repair",
                                           reason=f"{len(bad)} contact(s) would be read as covalent bonds by Meeko")
                else:
                    final.write_text(prot_text)
                text, n = _strip_clashing_hydrogens(final.read_text(), rep)
                if n:
                    final.write_text(text)
                    rep.warn(f"{n} hydrogen(s) removed because they sit on top of other atoms and could not be "
                             "repositioned without moving heavy atoms (listed in the Prep Report)")
        _clashes(final.read_text(), rep)
        for label in stagecheck.missing_amide_hydrogens(final.read_text()):
            rep.warn(f"{label} has no backbone hydrogen (chain start after a gap): its N-H donor is missing")
        lacking = stagecheck.missing_sidechain_hydrogens(final.read_text())
        if lacking:
            rep.warn(f"{len(lacking)} residue(s) lack polar side-chain hydrogens that PDB2PQR did not place "
                     f"(donors missing from the receptor), first: {', '.join(lacking[:6])}")

        # merge kept hetero groups back for Meeko
        if het_out:
            het = gemmi.Structure()
            hm = gemmi.Model("1")
            for cname in {c for c, _ in het_out}:
                hc = gemmi.Chain(cname)
                for c2, r in het_out:
                    if c2 == cname:
                        hc.add_residue(r)
                hm.add_chain(hc)
            het.add_model(hm)
            lines = [l for l in final.read_text().splitlines() if l.startswith(("ATOM", "TER"))]
            lines += [l for l in het.make_pdb_string().splitlines() if l.startswith("HETATM")]
            (td / "for_meeko.pdb").write_text("\n".join(lines) + "\nEND\n")
            meeko_in = td / "for_meeko.pdb"
        else:
            meeko_in = final
        cancel.check()
        defer_het = opts.pdbqt_writer != "meeko" and not opts.flex_residues
        meeko_in = _drop_untemplated(meeko_in, td / "meeko_filtered.pdb", rep, defer=defer_het)
        meeko_in = _drop_malformed(meeko_in, td / "meeko_sane.pdb", rep)
        rigid, flex = _write_pdbqt(meeko_in, outdir / "receptor", opts, rep)
        stagecheck.check_pdbqt(meeko_in.read_text(), rigid, rep.dropped, flex)

    _annotate_dropped(rep, model, site_pos)
    final_model = gemmi.read_pdb_string(final.read_text())[0]
    atoms_after = sum(len(r) for ch in final_model for r in ch)
    charge = 0
    for ch in final_model:
        for r in ch:
            charge += CHARGE.get(r.name, 0)
    # Termini (+1/-1 per chain) cancel, so the estimate is residue-based only.
    report = {
        "source": src.name, "options": {k: v for k, v in vars(opts).items()},
        "changes": rep.changes, "warnings": rep.warnings, "tools": rep.tools,
        "warning_details": [{"message": w, "severity": warning_severity(w)} for w in rep.warnings],
        "atoms_before": atoms_before, "atoms_after": atoms_after,
        "net_charge_estimate": charge, "chains": chains, "list_version": hetatm.LIST_VERSION,
        "site_residues": sorted(site_res_ids),
        # "degraded": the PDBQT is missing residues. The app must show this prominently (principle 2).
        # Left-out residues far from the site (sugars, stubs) do not change what is docked: noted, status stays ok.
        "status": "degraded" if any(_is_degrading(d) for d in rep.dropped) else "ok",
        "dropped_residues": rep.dropped,
        "assessment": assessment, "pdbqt": rep.pdbqt_info, "chain_map": chain_map,
    }
    return ReceptorResult(final, rigid, flex, ligands, report)
