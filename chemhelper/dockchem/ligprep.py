"""Ligand preparation pipeline (plan §6).

Stages: sanitize -> standardize -> protonate -> tautomers -> stereo -> 3D -> Meeko -> checks.
Every transformation is appended to `steps` so nothing happens silently (principle 2).
"""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import AllChem, QED, Descriptors, Lipinski, rdMolDescriptors
from rdkit.Chem.FilterCatalog import FilterCatalog, FilterCatalogParams
from rdkit.Chem.MolStandardize import rdMolStandardize
from rdkit.Chem.EnumerateStereoisomers import EnumerateStereoisomers, StereoEnumerationOptions

from .errors import InvalidInput
from .ligand_templates import is_untrusted

# Torsions above this make the Vina search less reliable (plan §9.3 "high torsion count").
TORSION_WARN = 15
# Elements Vina/Meeko handle without special care.
COMMON_ELEMENTS = {"H", "C", "N", "O", "F", "P", "S", "Cl", "Br", "I"}

_PARAMS = FilterCatalogParams()
_PARAMS.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS)
_PARAMS.AddCatalog(FilterCatalogParams.FilterCatalogs.BRENK)
_ALERTS = FilterCatalog(_PARAMS)


@dataclass
class LigandPrepOptions:
    strip_salts: bool = True
    neutralize: bool = True
    protonate: bool = True
    ph: float = 7.4
    ph_range: float = 0.0            # +/- around ph; 0 = dominant state only
    tautomers: str = "keep"          # keep | canonical | enumerate
    stereo: str = "keep"             # keep | enumerate  (enumerate only undefined centres)
    max_variants: int = 8            # cap for each enumeration
    seed: int = 42
    n_conformers: int = 1
    keep_3d: bool = True             # preserve input coordinates (needed for redocking RMSD)
    rigid_macrocycles: bool = False
    pdbqt: bool = True


@dataclass
class PreparedLigand:
    name: str
    status: str                      # prepared | warning | failed
    smiles: str = ""
    inchikey: str = ""
    mol: Chem.Mol | None = None      # 3D, explicit Hs
    pdbqt: str = ""
    steps: list[str] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)
    error: str = ""
    properties: dict = field(default_factory=dict)

    def summary(self) -> dict:
        return {"name": self.name, "status": self.status, "smiles": self.smiles,
                "inchikey": self.inchikey, "steps": self.steps, "warnings": self.warnings,
                "error": self.error, "properties": self.properties}


def properties(mol: Chem.Mol) -> dict:
    """Property panel (plan §6.2). `mol` may carry Hs; counts use the heavy-atom graph."""
    m = Chem.RemoveHs(mol)
    mw, logp = Descriptors.MolWt(m), Descriptors.MolLogP(m)
    hbd, hba = Lipinski.NumHDonors(m), Lipinski.NumHAcceptors(m)
    tpsa, rotb = rdMolDescriptors.CalcTPSA(m), rdMolDescriptors.CalcNumRotatableBonds(m)
    lip_viol = sum([mw > 500, logp > 5, hbd > 5, hba > 10])
    alerts = [e.GetDescription() for e in _ALERTS.GetMatches(m)]
    return {
        "mw": round(mw, 2), "clogp": round(logp, 2), "hbd": hbd, "hba": hba,
        "tpsa": round(tpsa, 1), "rotatable_bonds": rotb,
        "heavy_atoms": m.GetNumHeavyAtoms(), "qed": round(QED.qed(m), 3),
        "lipinski_violations": lip_viol, "lipinski_pass": lip_viol <= 1,   # one violation allowed
        "veber_pass": rotb <= 10 and tpsa <= 140,
        "alerts": alerts,                       # flags only, never filters (plan §6.2)
    }


def _has_3d(mol: Chem.Mol) -> bool:
    if mol.GetNumConformers() == 0:
        return False
    conf = mol.GetConformer()
    return conf.Is3D() and bool((abs(conf.GetPositions()[:, 2]) > 1e-4).any())


def _standardize(mol, opts, steps):
    if opts.strip_salts and len(Chem.GetMolFrags(mol)) > 1:
        n = len(Chem.GetMolFrags(mol))
        mol = rdMolStandardize.LargestFragmentChooser(preferOrganic=True).choose(mol)
        steps.append(f"standardize: kept largest of {n} fragments")
    if opts.neutralize:
        neutral = rdMolStandardize.Uncharger().uncharge(mol)
        if Chem.MolToSmiles(neutral) != Chem.MolToSmiles(mol):
            steps.append("standardize: neutralized charges")
        mol = neutral
    return mol


def _protonate(mol, opts, steps, warnings):
    from dimorphite_dl import protonate_smiles
    smi = Chem.MolToSmiles(mol)
    lo, hi = opts.ph - opts.ph_range, opts.ph + opts.ph_range
    try:
        # precision=1.0 (one pKa std-dev) returns every state that is plausible at this pH. Tighter
        # settings over-ionise borderline groups (e.g. a para-acylamino phenol, pKa ~9.5).
        states = protonate_smiles(smi, ph_min=lo, ph_max=hi, precision=1.0,
                                  max_variants=opts.max_variants)
    except Exception as e:  # noqa: BLE001 - external tool
        warnings.append(f"protonation failed ({e!r}); kept input state")
        return [mol]
    out = []
    for s in states or []:
        m = Chem.MolFromSmiles(s)
        if m is not None:
            out.append(m)
    if not opts.ph_range and len(out) > 1:
        # Dominant-state mode: when ionisation is ambiguous, keep the input (neutral) state rather
        # than guess; the user can enumerate states with a pH range.
        if smi in {Chem.MolToSmiles(m) for m in out}:
            warnings.append("ionisation state ambiguous at pH %g; kept the input state "
                            "(alternatives: %s)" % (opts.ph, ", ".join(sorted(Chem.MolToSmiles(m) for m in out if Chem.MolToSmiles(m) != smi))))
            steps.append("protonation: ambiguous, kept input state (Dimorphite-DL)")
            return [mol]
    if not out:
        warnings.append("protonation returned no valid states; kept input state")
        return [mol]
    steps.append(f"protonation: {len(out)} state(s) at pH {lo:g}-{hi:g} (Dimorphite-DL)"
                 if opts.ph_range else f"protonation: dominant state at pH {opts.ph:g} (Dimorphite-DL)")
    return out


def _tautomers(mol, opts, steps):
    if opts.tautomers == "keep":
        return [mol]
    te = rdMolStandardize.TautomerEnumerator()
    if opts.tautomers == "canonical":
        steps.append("tautomer: canonical")
        return [te.Canonicalize(mol)]
    if opts.tautomers == "enumerate":
        res = list(te.Enumerate(mol))[: opts.max_variants] or [mol]
        steps.append(f"tautomer: enumerated {len(res)} (cap {opts.max_variants})")
        return res
    raise InvalidInput(f"unknown tautomer mode {opts.tautomers}")


def _stereo(mol, opts, steps, warnings):
    undefined = Chem.FindMolChiralCenters(mol, includeUnassigned=True, useLegacyImplementation=False)
    undefined = [i for i, c in undefined if c == "?"]
    if not undefined or opts.stereo == "keep":
        if undefined:
            warnings.append(f"{len(undefined)} undefined stereocentre(s) left as-is")
        return [mol]
    isos = list(EnumerateStereoisomers(
        mol, options=StereoEnumerationOptions(onlyUnassigned=True, maxIsomers=opts.max_variants,
                                              rand=opts.seed)))
    steps.append(f"stereo: enumerated {len(isos)} isomer(s) for {len(undefined)} undefined centre(s)")
    if len(isos) >= opts.max_variants:
        warnings.append(f"stereoisomer enumeration capped at {opts.max_variants}")
    return isos


def _embed(mol, opts, steps, warnings):
    """Return (mol with Hs and conformer(s), ok)."""
    mh = Chem.AddHs(mol, addCoords=_has_3d(mol))
    if _has_3d(mol) and opts.keep_3d:
        steps.append("3D: preserved input coordinates")
        return mh
    ps = AllChem.ETKDGv3()
    ps.randomSeed = opts.seed
    cids = list(AllChem.EmbedMultipleConfs(mh, opts.n_conformers, ps))
    if not cids:
        ps.useRandomCoords = True
        cids = list(AllChem.EmbedMultipleConfs(mh, opts.n_conformers, ps))
    if not cids:
        raise InvalidInput("3D embedding failed")
    if AllChem.MMFFHasAllMoleculeParams(mh):
        AllChem.MMFFOptimizeMoleculeConfs(mh, maxIters=2000)
        steps.append(f"3D: ETKDGv3 (seed {opts.seed}) + MMFF94, {len(cids)} conformer(s)")
    else:
        AllChem.UFFOptimizeMoleculeConfs(mh, maxIters=2000)
        warnings.append("MMFF94 parameters unavailable; used UFF")
        steps.append(f"3D: ETKDGv3 (seed {opts.seed}) + UFF, {len(cids)} conformer(s)")
    return mh


def _to_pdbqt(mol, opts, warnings):
    from meeko import MoleculePreparation, PDBQTWriterLegacy
    kw = {"rigid_macrocycles": True} if opts.rigid_macrocycles else {}
    setups = MoleculePreparation(**kw).prepare(mol)
    text, ok, err = PDBQTWriterLegacy.write_string(setups[0])
    if not ok:
        raise InvalidInput(f"Meeko failed: {err}")
    ntors = text.count("BRANCH") // 2
    if ntors >= TORSION_WARN:
        warnings.append(f"{ntors} torsions: search may be unreliable; consider higher exhaustiveness")
    return text


def prepare_ligand(source: str | Chem.Mol, name: str = "ligand",
                   opts: LigandPrepOptions | None = None) -> list[PreparedLigand]:
    """Prepare one input; may return several variants (states/tautomers/stereoisomers).
    Never raises for chemistry problems: failures come back as status='failed'."""
    opts = opts or LigandPrepOptions()
    steps: list[str] = []
    warnings: list[str] = []
    try:
        if isinstance(source, Chem.Mol) and is_untrusted(source):
            return [PreparedLigand(name, "failed", error=(
                "bond orders are unknown (the molecule came from a PDB-style file, which has none). Preparing "
                "it would invent chemistry. Establish bond orders with ligand.from_structure (SMILES, bundled "
                "dictionary, or net charge for files with hydrogens)."))]
        mol = Chem.MolFromSmiles(source) if isinstance(source, str) else Chem.Mol(source)
        if mol is None:
            return [PreparedLigand(name, "failed", error="RDKit could not parse/sanitize input")]
        if isinstance(source, Chem.Mol):
            Chem.SanitizeMol(mol)
        had_3d = _has_3d(mol) and opts.keep_3d
        mol = _standardize(mol, opts, steps)
        if opts.protonate and not had_3d:
            states = _protonate(mol, opts, steps, warnings)
        else:
            if opts.protonate and had_3d:
                steps.append("protonation: skipped to preserve input 3D coordinates")
            states = [mol]
        variants = []
        for s in states:
            for t in _tautomers(s, opts, steps):
                variants.extend(_stereo(t, opts, steps, warnings))
        results, seen = [], set()
        for i, v in enumerate(variants):
            smi = Chem.MolToSmiles(v)
            if smi in seen:
                continue
            seen.add(smi)
            vname = name if len(variants) == 1 else f"{name}_{i + 1}"
            w, st = list(warnings), list(steps)
            try:
                m3 = _embed(v, opts, st, w)
                elems = {a.GetSymbol() for a in m3.GetAtoms()} - COMMON_ELEMENTS
                if elems:
                    w.append(f"uncommon element(s) {sorted(elems)}: Vina typing may be unreliable")
                pq = _to_pdbqt(m3, opts, w) if opts.pdbqt else ""
                results.append(PreparedLigand(
                    vname, "warning" if w else "prepared", smi, Chem.MolToInchiKey(v), m3, pq,
                    st, w, properties=properties(m3)))
            except Exception as e:  # noqa: BLE001
                results.append(PreparedLigand(vname, "failed", smi, steps=st, warnings=w, error=str(e)))
        return results
    except Exception as e:  # noqa: BLE001
        return [PreparedLigand(name, "failed", steps=steps, warnings=warnings, error=str(e))]


def read_records(path: str | Path):
    """Yield (name, mol_or_smiles_or_None, error) from an SDF or SMILES file."""
    path = Path(path)
    ext = path.suffix.lower()
    if ext in (".sdf", ".mol"):
        for i, m in enumerate(Chem.SDMolSupplier(str(path), removeHs=False, sanitize=True)):
            if m is None:
                yield f"record_{i + 1}", None, "unreadable SDF record"
            else:
                yield (m.GetProp("_Name") if m.HasProp("_Name") and m.GetProp("_Name") else f"record_{i + 1}"), m, ""
    elif ext in (".smi", ".smiles", ".txt"):
        for i, line in enumerate(path.read_text().splitlines()):
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.split(None, 1)
            yield (parts[1].strip() if len(parts) > 1 else f"mol_{i + 1}"), parts[0], ""
    else:
        raise InvalidInput(f"unsupported library format {ext}")


def prepare_library(records, opts: LigandPrepOptions | None = None) -> dict:
    """Library mode (plan §6.3): dedupe by InChIKey, per-molecule status, failures never vanish."""
    rows, prepared, seen = [], [], {}
    for name, src, err in records:
        if src is None:
            rows.append({"name": name, "status": "failed", "error": err, "warnings": [], "steps": []})
            continue
        probe = Chem.MolFromSmiles(src) if isinstance(src, str) else src
        if probe is not None:
            key = Chem.MolToInchiKey(Chem.RemoveHs(probe))
            if key in seen:
                rows.append({"name": name, "status": "duplicate", "error": f"duplicate of {seen[key]}",
                             "warnings": [], "steps": []})
                continue
            seen[key] = name
        for p in prepare_ligand(src, name, opts):
            rows.append(p.summary())
            if p.status != "failed":
                prepared.append(p)
    counts = {s: sum(r["status"] == s for r in rows) for s in ("prepared", "warning", "failed", "duplicate")}
    return {"rows": rows, "counts": counts, "prepared": prepared}
