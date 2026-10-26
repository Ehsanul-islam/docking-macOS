# 1000-structure survey: what went wrong, where and why

Run date 2026-10-25/26. Software at commit `6f6d1d7` (nothing changed during the run; only run tooling was touched).
Set: 1002 entries, the most-studied proteins (largest UniProt groups) per organism group, best-resolution representative,
X-ray <= 2.5 A / cryo-EM <= 3.5 A / NMR; <= 100,000 atoms. Manifest `tools/batch/sets/most_studied_1000.json`.
Each entry: assessment, then full receptor preparation with its largest drug-like ligand as the site (333 entries have
no ligand, so no site). Per-entry limit 15 min, 4 parallel workers. **No fixes were made; this is a problem catalogue.**
Raw per-entry data: `tools/batch/survey_notes.py` over the run's `proteins.jsonl` (not committed; regenerate with the batch tool).

## Headline

| Outcome | Entries | Meaning |
|---|---|---|
| ok | 864 (86.2%) | prepared, nothing material left out (or only far from the site) |
| degraded | 79 (7.9%) | prepared with something left out, or no site to judge it against |
| unusable (refused with reason) | 36 (3.6%) | 32 metal cofactor at the site, 4 flagged distorted (**all 4 are false positives**, P4) |
| timeout (15 min) | 14 (1.4%) | P5 |
| failed | 9 (0.9%) | 8 non-standard residues (P3), 1 chain-name limit (P8) |

Prepared (ok + degraded): 943 / 1002 = 94.1%. Built-in writer produced the PDBQT in 941; Meeko in 2. Meeko cross-check:
922 agree, 4 disagree, 15 **Meeko itself failed** (the table writer rescued these), 2 n/a.
By group (ok/degraded/unusable/failed/timeout): human 370/13/5/4 (+ timeouts), animal 124/17/6/2, bacteria 206/26/16/0,
plant 36/10/3/0, other 126/13/6/3. By method: X-ray 781 ok of 896; cryo-EM 52 ok of 60 (+11 of the 14 timeouts are cryo-EM);
NMR 29 of 30. Bacteria and plants refuse/degrade most (metalloproteins, chlorophyll, Fe-S).

## Problems, ordered by how many entries they affect

### P1. Metal cofactors at the site: 32 refused (3.2%), plus 32 HEM residues in degraded entries
- **What:** HEM 13, 4Fe-4S (SFN) 4, chlorophyll/bacteriochlorophyll (CLA, BCL) 4, Co 3, Cu 4, Ni 2, W 1, Mo 1.
- **Where:** site-aware check at the start of `prepare_receptor` (`_unsupported_reason`). **Why:** no verified AutoDock
  type for porphyrins/clusters, and Co/Cu/Ni/W/Mo/Cd/Hg have no type in Meeko's AD4 table.
- **Note:** this is the designed behaviour, it is also the largest single class of refusals. Bacteria 16, animal 6, other 6.
- **Candidate fix (later):** type common metal cofactors (needs validated types); confirm whether Vina itself has types for
  Cu/Co/Ni (Meeko's table does not). Decision for the owner.

### P2. "Degraded" is mostly "no site defined": 74 of 79 degraded entries have no reference ligand
- **What:** with no site, nothing left out can be judged near or far, so any omission marks the entry degraded.
- **Left out in those entries:** NAG sugars 90 residues, Cu ions 43, HEM 32, Cd 12, Hg 11, Co 9, PLP 7, Ni 5, residues with no
  backbone 16. **Where:** `_annotate_dropped` returns `at_site=None` without a site; status rule treats None as degraded.
- **Why it matters:** blind docking (no ligand) is a normal use; those results look worse than they are, and a heme anywhere
  in the structure is not refused because there is no site to be near.
- **Candidate fix:** define "site" for blind docking (whole-protein or user box), judge omissions against it; show far/near
  in the report and UI.

### P3. Non-standard residues fail the hydrogen stage check: 8 entries hard-fail
- **Entries/residues:** MLY (dimethyl-lysine) 3F9X, 6BHD, 9D1L; UNK (unidentified) 3ED7, 9NVO, 7Z0H; SEC (selenocysteine) 7A0D;
  ASX (Asp/Asn ambiguous) 3E2O.
- **Where:** `stagecheck.check_hydrogens` after PDB2PQR ("X has no hydrogens"). **Why:** PDB2PQR keeps residues it has no
  template for but does not protonate them; my check treats any residue gemmi calls an amino acid as needing hydrogens.
- **Inconsistency:** the same kind of residue is *dropped by name* elsewhere (no Meeko template), but here it fails the whole
  entry. Expected behaviour is probably drop-by-name + degraded (or map MLY->LYS, SEC->CYS with a warning).
- **Note:** `assess` does not flag these in advance (`UNK` only if > 10% of residues).

### P4. Deuterated (neutron) structures are refused as "distorted": 4 entries, all false positives
- **Entries:** 9YFR, 6S2M, 5AE0, 5NW3 (sub-1.1 A resolution neutron/X-ray joint refinements).
- **Where/why:** `assess._distorted`: deuterium atoms (element `D`) are not recognised as hydrogen, so O-D at 0.84 A counts
  as two overlapping heavy atoms. Same defect class as the earlier O-H false positive (which was fixed for `H` only).
- **Check also:** other code paths that remove hydrogens (`st.remove_hydrogens`, PDB2PQR input) were not exercised on D;
  untested.

### P5. Timeouts: 14 entries (1.4%), all large (1,350 to 12,200 residues) and 13 of 14 without a site
- **Entries:** cryo-EM 10QQ, 9BZ0, 9BLZ, 7NVL, 9QGG, 8SH9, 7Z8R, 8BKZ, 7NK7, 6RDB, 6RD4, 9NPX; X-ray 1JZ7, 3SXX.
- **Profiled one (9QGG, 10.7k atoms, 1,350 residues), run alone:** assessment 1 s, PDBFixer 1 s, PDB2PQR 20 s, close-contact
  detection 2 s, **hydrogen-repair minimisation (OpenMM) 490 s**, Meeko cross-check 37 s; total 555 s. In the batch it exceeded
  900 s: four concurrent jobs share four cores (likely OpenMM thread oversubscription; not separately measured).
- **Where:** `recprep._minimize` (restrained minimisation used to repair hydrogen clashes, triggered whenever clashes exist).
  **Why:** cost grows faster than linearly with residue count; no timeout, progress or size limit inside the helper.
- **Also:** p99 time 389 s; median 12 s. The Meeko cross-check adds 37 s here and more on bigger entries.
- **Candidate fixes:** minimise only the clashing neighbourhood; limit threads per job; skip/shorten the cross-check above a size;
  per-stage timeouts with a clear message and `degraded`/`failed` instead of an opaque timeout.

### P6. Meeko fails on copper-containing structures; the table writer then drops the copper
- **What:** 15 entries where the Meeko cross-check failed ("atom_type is None"; checked 5 of 6 sampled: all contain Cu).
  Meeko has a CU residue template but no AutoDock type for copper.
- **Effect on us:** those 15 are the cases where the built-in writer is demonstrably better than the old route. But the table
  writer has no rule for CU either, so 43 Cu residues are left out ("no atom-type rules for CU") and the entries degrade.
- **Candidate fix:** decide whether Vina can score Cu at all (needs the Vina type table, not Meeko's), then add a rule or
  refuse clearly.

### P7. Table writer incomplete for rare residue states: 2 entries fall back to Meeko
- 1PJX (LYS43), 3F9M (GLU4): "no rules for". Likely atom naming or protonation variant not in the harvested table (corpus
  came from the earlier ~46 preparations). When Meeko takes over, cofactors Meeko cannot template (3F9M: GLC) are left out
  (the dictionary route only runs on the table path).
- **Candidate fix:** harvest the table from this larger corpus; apply the dictionary route to Meeko-written outputs too.

### P8. Chain-name aliasing runs out: 1 entry (8TOX, cryo-EM, many multi-letter chain IDs)
- `_alias_chains` needs a free one-letter name per multi-letter chain (62 total). **Fix idea:** ask the user to pick chains,
  or write mmCIF end to end (the PDB format is the limit; also the 99,999-atom limit).

### P9. Writer vs Meeko disagreements: 4 entries (6R2W, 4QYN, 8UGP, 8RUP)
- Pattern: one residue where Meeko writes more polar hydrogens (HD) than the table path: TYR73/TYR264 (HD 1 vs 2), ARG134
  (5 vs 6), CYX17 (1 vs 3). Looks like residues at chain starts after a gap (no amide/ammonium hydrogens from PDB2PQR; the
  existing warning "has no backbone hydrogen" covers the amide case). Not diagnosed per residue.

### P10. Frequent warnings that may be noise or may hide a science issue
- "unusual protonation" for HIP (230+ entries), GLH, ASH, LYN, CYM; metals retained: Zn 59, Mg 57.
- **Cysteine:** 59 entries get a deprotonated Cys (CYM, PROPKA) represented as CYX in the PDBQT, so "its negative charge is
  not modelled" (59 matching warnings). Whether PROPKA's thiolates are plausible at pH 7.4 for these entries was not checked.
- The volume of warnings makes them hard to act on; a severity split is probably needed (hold for the UI work).

### P11. Pre-processing volume (information only)
Alternate conformations in 762 assessments, missing atoms in 364, backbone stubs 35, distorted residues 8 (real). These are
handled; they show how rarely a PDB entry is "ready" as downloaded (about 12% here).

## What worked
- No crashes of the helper: every entry ended in a classified outcome (apart from the two harness bugs below).
- The built-in writer produced 941 PDBQTs; Meeko could not prepare 15 of them at all.
- Dictionary-typed cofactors (FAD, NAP, FMN, ...) and ions went through without incident.

## Caveats of this run
- Two entries (3F9M, 1PJX) first crashed in the run tooling (a string where a dict was expected); fixed in the tool, re-run.
- The container was reclaimed twice (session idle); the runner resumed from its results file, so no entry was run twice.
- Timeouts are 15 min with 4 jobs sharing 4 cores; some "timeouts" may finish alone (9QGG did, in 555 s).
- Reference ligand = the largest drug-like group by the app's own classification, which can pick a cofactor or a buffer-like
  ligand; docking accuracy was **not** measured on these entries (preparation only).
- The set is skewed toward well-studied proteins and excludes entries over 100,000 atoms, so very large complexes are
  under-represented.
