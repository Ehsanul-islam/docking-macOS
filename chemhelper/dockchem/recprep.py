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

from . import hetatm
from .errors import InvalidInput, ToolFailure
from .provenance import record

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
    minimize: bool = False
    flex_residues: list[str] = field(default_factory=list)   # e.g. ["A:15"]
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

    def change(self, stage, residue, change, reason, tool):
        self.changes.append({"stage": stage, "residue": residue, "change": change,
                             "reason": reason, "tool": tool})

    def warn(self, msg):
        if msg not in self.warnings:
            self.warnings.append(msg)


def _rid(chain, res) -> str:
    return f"{chain}:{res.name}{res.seqid.num}{res.seqid.icode.strip()}"


def _heavy(res) -> int:
    return sum(1 for a in res if a.element.name != "H")


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


def _run(cmd, cwd):
    import os
    return subprocess.run(cmd, cwd=cwd, capture_output=True, text=True, env={**os.environ, **OFFLINE_ENV})


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
    worst = max((math.dist(b[k], a[k]) for k in b if k in a), default=0.0)
    if worst > MAX_HEAVY_SHIFT:
        return f"heavy atom moved {worst:.2f} A (limit {MAX_HEAVY_SHIFT})"
    for (ch, num), _ in {(k[0], k[1]): 0 for k in a}.items():
        for (n1, n2), (lo, hi) in BOND_RANGE.items():
            p, q = a.get((ch, num, n1)), a.get((ch, num, n2))
            if p and q and not lo <= math.dist(p, q) <= hi:
                return f"implausible {n1}-{n2} bond ({math.dist(p, q):.2f} A) at {ch}:{num}"
    return ""


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
                             mm.Platform.getPlatformByName("CPU"))
        sim.context.setPositions(pdb.positions)
        sim.minimizeEnergy(maxIterations=500)
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


def _drop_untemplated(src: Path, dst: Path, rep) -> Path:
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
                continue
        lines.append(line)
    for ch, num in sorted(renamed):
        rep.warn(f"{ch}:CYS{num} has no HG (thiolate or unbridged): represented as CYX in the PDBQT, "
                 "so its negative charge is not modelled")
    for (ch, num, name), n in dropped.items():
        rep.warn(f"{ch}:{name}{num} left out of the PDBQT: no bundled Meeko template "
                 f"({n} atoms); Meeko would otherwise download it from RCSB")
        rep.dropped.append({"residue": f"{ch}:{name}{num}", "reason": "no bundled Meeko template"})
    if not dropped and not renamed:
        return src
    dst.write_text("\n".join(lines) + "\n")
    return dst


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
    if r.returncode != 0:
        keys = set(failing(r))
        lines, stripped = [], set()
        for l in pdb_path.read_text().splitlines():
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
    n_h = sum(1 for ch in model for r in ch for a in r if a.element.name == "H")
    atoms_before = sum(len(r) for ch in model for r in ch)
    if n_h:
        st.remove_hydrogens()
        rep.change("hydrogens", "all", f"removed {n_h} input hydrogens", "re-added consistently by PDB2PQR", "gemmi")
    model = st[0]

    # site residue occupancy / B-factor diagnostics need pre-collapse info
    partial_occ = {_rid(ch.name, r) for ch in model for r in ch if any(a.occ < 1.0 for a in r)}
    _collapse_altlocs(model, rep)

    refs_all = _find_reference(model, opts.reference_ligand)
    refs = refs_all[:1]
    if len(refs_all) > 1:
        rep.warn(f"{len(refs_all)} copies of {opts.reference_ligand} found; site defined by "
                 f"{_rid(refs[0][0], refs[0][1])}; others: "
                 + ", ".join(_rid(c, r) for c, r in refs_all[1:]) + " (plan §8.2: pick one)")
    ref_pos = [a.pos for _, r in refs for a in r]
    if opts.reference_ligand and not refs:
        raise InvalidInput(f"reference ligand {opts.reference_ligand} not found")

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

    # build polymer-only structure for PDBFixer / PDB2PQR
    poly = gemmi.Structure()
    poly.spacegroup_hm, poly.cell, poly.name = st.spacegroup_hm, st.cell, st.name
    pm = gemmi.Model("1")
    for ch in model:
        if ch.name not in chains:
            continue
        pc = gemmi.Chain(ch.name)
        for res in ch:
            if _is_polymer(res):
                if res.name in hetatm.MODIFIED and opts.nonstandard == "remove":
                    rep.change("nonstandard", _rid(ch.name, res), "removed", "user chose to remove", "gemmi")
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
        fixed = _pdbfixer(poly.make_pdb_string(), opts, rep)
        (td / "fixed.pdb").write_text(fixed)
        rep.tools.append(record("pdbfixer", _ver("pdbfixer"), {"rebuild_missing_atoms": opts.rebuild_missing_atoms}))
        before = _residue_names(fixed)
        prot = _pdb2pqr(td / "fixed.pdb", opts, rep, td)
        prot_text = prot.read_text()
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
                rep.change("protonation", f"{key[0]}:{old}{key[1]}", f"{old} -> {new}",
                           f"PROPKA pKa at pH {opts.ph:g}" if new != "HID" and new != "HIE" else
                           "histidine tautomer from H-bond network", "PDB2PQR+PROPKA")
                if new in UNUSUAL_STATES:
                    rep.warn(f"unusual protonation {key[0]}:{new}{key[1]} ({UNUSUAL_STATES[new]})")
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
        meeko_in = _drop_untemplated(meeko_in, td / "meeko_filtered.pdb", rep)
        rigid, flex = _meeko(meeko_in, outdir / "receptor", opts, rep)

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
        "atoms_before": atoms_before, "atoms_after": atoms_after,
        "net_charge_estimate": charge, "chains": chains, "list_version": hetatm.LIST_VERSION,
        "site_residues": sorted(site_res_ids),
        # "degraded": the PDBQT is missing residues. The app must show this prominently (principle 2).
        "status": "degraded" if rep.dropped else "ok", "dropped_residues": rep.dropped,
    }
    return ReceptorResult(final, rigid, flex, ligands, report)
