# Phase 2 status (protein & ligand preparation)

Exit criteria (plan §18): §12.2 reference set passes; round-trip tests pass.

| Plan item | State |
|---|---|
| §6 ligand pipeline (sanitize, salts, neutralise, pH, tautomers, stereo, ETKDGv3+MMFF, Meeko, checks) | Done, tested |
| §6.2 property panel (MW, cLogP, HBD/HBA, TPSA, RotB, QED, Lipinski, Veber, PAINS/Brenk flags) | Done |
| §6.3 library mode (InChIKey dedupe, per-molecule status, failures listed) | Done |
| §5 receptor pipeline (model/chain/altloc/water/HETATM, PDBFixer, PDB2PQR+PROPKA, flips, OpenMM restrained min, Meeko, flex residues, AlphaFold pLDDT) | Done, tested on synthetic structures |
| §5.2 Prep Report (change table, warnings, atom counts, tool versions) | Done (`prep_report.json`) |
| §7 format matrix (PDB, mmCIF, PDBQT, SDF/MOL, MOL2, SMILES, XYZ), lossy warnings, pose PDBQT->SDF | Done, round-trip tested |
| RPC: `ligand.prepare`, `library.prepare`, `receptor.prepare`, `format.convert`, `format.sniff`, `poses.to_sdf` | Done |

## Not done / needs a decision

- **§12.2 reference set (50-100 expert-reviewed PDB entries) is NOT done.** RCSB is unreachable from the dev
  sandbox, so everything is tested on synthetic peptides. Real-structure validation is required before the
  Phase 2 exit criterion can be claimed.
- PQR import/export, per-residue manual His/Asp/Glu/Lys/Cys state overrides, per-flip accept/reject: not implemented.
- Full CCD bundle for HETATM classification: curated lists + heavy-atom count fallback only (`hetatm.py`).
- Kept cofactors/metals/glycans get no hydrogens and are dropped from PDBQT by Meeko if it lacks a template (reported).
- D-3 (flip tool) provisionally resolved as PDB2PQR's optimiser; Reduce not evaluated. Licences (PROPKA, Meeko, OpenMM, PDB2PQR) still [VERIFY].
- MOL2 writer is our own (RDKit has none): Gasteiger charges, Sybyl types for common elements only.

## Findings worth knowing

- Dimorphite-DL at tight precision over-ionises borderline groups (paracetamol phenol). Default mode now keeps the
  input state and warns when ionisation is ambiguous; a pH range enumerates states.
- If PDB2PQR cannot debump a clash, Meeko fails with a valence error. The report surfaces the clash warning and the
  Meeko error carries a hint.
- Default box/prep defaults here are unbenchmarked (plan §12.4: needs advisor sign-off).

## Gap 3: bond orders for ligands taken from structure files (done on this branch)

**Bug found by audit:** ligands extracted from a PDB came out as all-single-bond heavy-atom skeletons and the
ligand pipeline prepared them anyway (imatinib 1IEP became a saturated molecule; biotin 1STP's acid became a
gem-diol), returning only a "warning". Fixed:

- Molecules read from PDB-style files are flagged *untrusted*. `prepare_ligand` and the SDF/MOL/MOL2/SMILES writers
  refuse them with a typed `NeedsTemplate` error (`needs`: `smiles` or `none`).
- `ligand.from_structure` (module `ligand_templates`) establishes bond orders from exactly one source:
  user SMILES > wwPDB Chemical Component Dictionary > perception from net charge (only if the file has hydrogens).
  Otherwise it raises `NeedsTemplate` saying what to supply. It never guesses.
- Template matching checks heavy-atom composition and connectivity, takes stereochemistry from the crystal
  geometry (carbon centres only; 3D "stereo" on phosphates is an artefact), adds hydrogens with coordinates,
  and audits double-bond lengths.
- `receptor.prepare` reports `bond_orders: resolved | needs_template` for each extracted ligand.
- **CCD bundle:** `data/ccd_ligands.json.gz` (51,493 components, 4 MB, flat SMILES + formula + name; built by
  `tools/ligand_dictionary/build_ccd.py` from components.cif.gz dated 2026-10-17). See docs/DATA_SOURCES.md.
- The earlier hand-written 13-entry dictionary was cross-checked against the CCD (all 13 identical after
  assignment) and then **removed**: one authoritative source is safer than two.

**Results on the 50 supplied structures** (`ligand_from_structure` over every non-additive, non-metal HETATM group):
27 of 62 ligand-in-entry cases resolve, all with a clean bond-length audit (imatinib, erlotinib, biotin, ATP, GTP,
NAD(P), cholesterol, 0WN, 0WM, F09, P0G, TXL, TA1 ...). The 35 refusals are deliberate:

| Reason | Examples | Why refused |
|---|---|---|
| Metal complex | HEM x9 | Fe coordination: bond orders ambiguous, Vina models metals poorly |
| Atoms missing from the file | PLM, 12P, DGA, RET, PLP, CRO | capping a truncated ligand would invent a different molecule |
| Covalently attached sugar (linking O absent) | NAG, MAN, BMA, GAL, FUC, FUL, NDG | not a free ligand; template has the extra glycosidic O |

**Limits, stated plainly**

- Only ligands whose defining atoms are all present are resolved. An explicit "dock this truncated ligand anyway"
  opt-in does not exist yet.
- The bond-length audit catches gross misassignment, not tautomer errors (C=O vs C=N are within each other's range).
- The CCD stores one reference protonation state per component; ionisation at a chosen pH is not applied to
  ligands that keep their crystal coordinates.
- The CCD snapshot is a point-in-time copy; the build script must be re-run to refresh it (record the date).
- A net charge alone cannot establish bond orders for a hydrogen-free ligand; the charge route only works when the
  file already has hydrogens.

## Gap 1: free selection and deletion (done on this branch)

- `dockchem/selection.py`: a selection language (chain, residue name/range, element, atom name, classes, B-factor,
  occupancy, `within`/`around`/`byres`, boolean logic) with located parse errors; see docs/SELECTION.md.
- Non-destructive `EditSession`: the command list is the state, so undo/redo/replay are exact and the input is never
  modified. `preview` reports atoms/residues affected plus warnings (partial residues, new chain breaks, pocket
  residues, majority/all-protein removal).
- RPC: `selection.validate`, `selection.preview`, `selection.apply` (stateless).
- Receptor prep: new `edits` option; each edit is recorded in the Prep Report. Edits that remove all protein or the
  reference ligand are refused with a clear message.
- Tests: 50, using real structures with independent oracles (brute-force distances, direct gemmi counts) and
  mutation-checked (breaking `around` or the partial-residue warning makes them fail).

Not done: the 3D highlight of a previewed selection and the linked sequence strip (viewer work, needs the Swift app).
