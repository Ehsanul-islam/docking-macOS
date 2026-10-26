# Reliability: how we find and fix problems at scale

Testing by hand on a few dozen structures does not converge: each structure exposed a new failure. The approach now:

1. **Remove fragile steps** (receptor PDBQT is written from a residue/atom-name type table; Meeko is a cross-check).
2. **Check at the door** (`structure.assess`: ready / fixable / unusable, with named reasons) and **after every stage**
   (`stagecheck.py`).
3. **Test at scale, automatically**, with resumable runners that survive hangs and crashes (`tools/batch/`), and turn
   every failure into a category, then a fix, then a regression test.

## What has been run (this branch)

| Corpus | Result |
|---|---|
| All 51,493 wwPDB Chemical Component Dictionary ligands through ligand prep + output audit | 50,000 of 50,009 in-scope ligands pass (99.98%). 1,115 contain metals/metalloids Vina cannot type and are refused with a message; 246 dictionary entries RDKit cannot parse; 123 single atoms/ions. |
| 300 widely used structures (Astex Diverse Set + evenly spaced PoseBusters entries), full receptor prep with the entry's ligand as site | 300/300 prepare: 270 ok, 30 degraded (residues named). Built-in PDBQT writer used for all; Meeko agrees on 297 (3 are residues that start after a chain gap, where PDB2PQR places no amide hydrogen and Meeko adds one; now warned by name). Median 9 s per entry, p95 69 s. |

Reproduce: `tools/batch/ligand_batch.py`, `tools/batch/protein_batch.py` (see their docstrings). The 300-entry list is
`tools/batch/sets/widely_used_300.json`.

## Bugs these runs found (all fixed, each with a regression test)

Ligands: hypervalent atoms break Gasteiger charges (zero charges + warning); force-field minimisation collapsing SF5
groups; torn 3D embeddings (9 A "N-S bond"); metals silently passing with unreliable types.

Proteins: assessment counted X-ray O-H at 0.84 A as distortion (49 of 300 wrongly refused); stage check rejected
PDB2PQR's normal equivalent-atom relabelling; His ring flips undetected; backbone-stub residues crashing PDB2PQR;
multi-letter mmCIF chain names; PDB2PQR crash on an N-terminal Pro without its ring; cofactors and sugars kept in
the site silently untyped.

## Known limits (stated, not hidden)

- 30 of 300 are "degraded": covalently attached glycans (NAG), metal cofactors (heme, Fe-S, single Co/Fe ions), a few
  residues with no backbone or impossible coordinates. Each left-out residue is listed by name with the reason.
- Hydrogens on kept cofactors come from the dictionary protonation state, not from the structure or pH.
- Partial charges are not written by the table writer (Vina ignores them); AutoDock4 scoring would need Meeko.
- Flexible side chains still require Meeko (they need a torsion tree).
- The corpus so far is docking benchmarks, which are cleaner than the PDB at large; the next corpora should include
  cryo-EM, NMR, large complexes, and old low-resolution entries.
- Structures judged "unusable" (Calpha-only, nucleic-acid-only, heavily distorted) are refused unless
  `allow_unusable` is set.
- Accuracy has been measured on one benchmark only (below); more is needed before any claim beyond it.

## Redocking accuracy (Astex Diverse Set, 85 complexes)

Method (`tools/benchmark/redock.py`): receptor prepared by the full pipeline from the structure downloaded from the
PDB; ligand starts from the benchmark's generated conformer (not the crystal pose) and goes through our ligand
preparation; box centred on the crystal ligand, extent + 6 A each side; unmodified AutoDock Vina 1.2.7 (official Linux
binary), exhaustiveness 8, fixed seed, 9 modes; success = top-ranked pose within 2.0 A heavy-atom RMSD (symmetry-aware).

| | Full pipeline | Standard route (benchmark's cleaned protein straight through Meeko) |
|---|---|---|
| Entries that could be prepared and docked | **85 / 85** | 32 / 85 (53 stop on unusual residues when offline) |
| Top-1 within 2 A | **49 / 85 = 57.6%** | 19 / 32 = 59.4% (of all 85: 22.4%) |
| Best of top 3 | 61 / 85 = 71.8% | 22 / 32 |
| Best of all 9 poses | 65 / 85 = 76.5% | 23 / 32 |
| Entries where both routes ran (paired) | top-1 19/32 = 59% | top-1 19/32 = 59% (best of 9: 24 vs 23) |

Reading: on identical inputs the pipeline matches the standard route, and it prepares the 53 complexes the standard
route cannot. The absolute rate (about 58% top-1) is what unmodified Vina typically achieves on this set; the
pipeline does not improve Vina's scoring and does not claim to. Receptors with something left out (12 entries, e.g.
heme or a covalent glycan): 6/12 top-1 versus 43/73 for clean receptors.

Two harness bugs found on the way (both fixed, which is why the first numbers looked worse): the receptor was built
around the wrong chain copy of the ligand; RDKit refuses to strip stereo-defining hydrogens in the benchmark's start
conformers, so four ligands could not be scored.

Remaining large misses (> 8 A): 1HVY and 1W1P (clean receptors, not yet diagnosed), 1MMV and 1R9O (heme left out).
Not yet run: the PoseBusters subset (215 entries in the 300-entry set), flexible-residue docking, other exhaustiveness.

## Helper hardening (cancellation, limits, fuzzing)

- **Cancellation works mid-job.** The server reads requests on one thread and runs them in order on a worker, so `cancel`
  is seen while a job runs; jobs check a token at stage boundaries, every 25 minimisation iterations and between library
  molecules, and helper subprocesses (PDB2PQR, Meeko) are killed (process group) on cancel. Before this, only the
  `test.sleep` demo could be cancelled.
- **Time limits:** each helper subprocess is stopped after 30 min (`DOCKCHEM_STAGE_TIMEOUT`); minimisation 4 min.
- **Input guards** (`limits.py`): missing/empty/folder/unreadable/over 300 MB (`DOCKCHEM_MAX_INPUT_MB`) files, stray bytes
  in coordinate records, wrong parameter types, ligands over 200 heavy atoms: all typed `InvalidInput`, never a crash.
- **`InternalError`** (-32099, with traceback) now means "a bug", distinct from expected rejections.
- **Fuzzing** (`tests/test_robustness.py`: 19 malformed-file kinds x pdb/cif, 60 seeded random mutations, malformed
  SMILES/selections/parameter types): found and fixed 5 bugs (sniff on a model-less file, stray bytes crashing assess,
  non-text parameters, lone-surrogate SMILES, and ligands with thousands of atoms hanging the helper).
