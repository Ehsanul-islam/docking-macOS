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
  refuse them with a typed `NeedsTemplate` error.
- `ligand.from_structure` (module `ligand_templates`) establishes bond orders from exactly one source:
  user SMILES > bundled dictionary > perception from net charge (only if the file already has hydrogens).
  Otherwise it raises `NeedsTemplate` saying what to supply. It never guesses.
- Template matching checks heavy-atom composition and connectivity, takes stereochemistry from the crystal
  geometry (carbon centres only; 3D "stereo" on phosphates is an artefact), adds hydrogens with coordinates,
  and audits double-bond lengths.
- `receptor.prepare` reports `bond_orders: resolved | needs_template` for each extracted ligand.
- Bundled dictionary: 13 entries (STI, AQ4, BTN, BEN, CAU, ATP, GTP, GDP, NAD, NDP, AP5, CLR, BU1). Each was
  validated against a real crystal ligand (`validated_on`); `tools/ligand_dictionary/build_dictionary.py` refuses to
  ship an entry that does not match. CI re-validates every entry.

**Limits, stated plainly**

- On your 50 entries 15 of 38 ligand cases resolve, but the dictionary was written from this same set, so that is
  *not* a coverage estimate for arbitrary PDB ligands. It will be much lower. A full CCD snapshot (public domain,
  wwPDB) is still needed; it could not be downloaded in the sandbox (RCSB blocked).
- Not supported, and reported as `needs_template`: metal complexes (HEM), covalently bound or incomplete ligands
  (RET, PLP, CRO), and anything not in the dictionary and without a user SMILES.
- A net charge alone cannot establish bond orders for a hydrogen-free ligand (there is nothing to perceive from);
  the earlier suggestion of a "net charge" fallback only works when the file has hydrogens.
- The bond-length audit catches gross misassignment, not tautomer errors (C=O vs C=N are within each other's range).
- The dictionary stores one reference protonation state per ligand; ionisation at a chosen pH is not applied to
  ligands that keep their crystal coordinates.
- Dictionary SMILES were written by hand and validated against structures, not against the CCD itself.
