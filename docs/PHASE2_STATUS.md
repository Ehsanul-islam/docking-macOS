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
