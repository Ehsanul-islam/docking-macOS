# 2D interaction images and one-panel quick prep

## 2D interaction image (`complex.image2d`, `dockchem/interactions.py`, `dockchem/render2d.py`)
Reproduces the figure pipeline supplied in `2D_interaction_images_bundle.zip` (RDKit ligand depiction, residue discs placed
by ray-casting then relaxed, coloured dashed lines by interaction type, distance labels, shared legend, matplotlib).
The profiler is our own (no PLIP/Open Babel): H-bonds (D...A <= 3.5 A, D-H...A >= 120 deg; rotatable OH/NH3+ donors
X-D...A >= 90 deg), salt bridges (<= 4.0 A), water bridges, hydrophobic (0.5-4.0 A), pi-stacking (<= 5.5 A, parallel or
T-shaped, offset <= 2.0 A), pi-cation (<= 6.0 A), van der Waals (<= 4.0 A). Every label comes from the profile.
Styles: `reference` (the supplied look), `mono` (print-friendly greys, dash patterns tell the types apart), `dark`.
Input: the prepared receptor PDB (it must have hydrogens) and a pose (PDBQT with Vina scores, or SDF), any pose index.
Example: `docs/examples/3ERT_OHT_docked_2D.png` (4-hydroxytamoxifen docked into ERalpha, Vina -8.7 kcal/mol);
`3ERT_OHT_native_2D_reference_style.png` (crystal pose: phenol OH to Glu353/Arg394, the known anchor).
Not verified: pixel-identity with the supplied images (font and layout come from the same code, but no common input exists);
the profiler's counts differ from PLIP's where the criteria are interpreted differently.

## Quick prep (`quickprep.run`, `dockchem/quickprep.py`)
One call: receptor file + ligand (SMILES, or a SDF/MOL/SMI/PDB file, one or many) -> assess, choose the site, prepare the
receptor, prepare the ligand(s), box, engine config. Presets: fast (exhaustiveness 8), standard (32, default), thorough
(64, 20 modes). Site: automatic (largest drug-like ligand), a named ligand, a box, or residues; no site = blind docking, said
plainly with a warning. Refuses unusable structures with the reasons. Returns every decision plus `vina.conf`, receptor
and ligand PDBQT files. Measured on 3ERT: 9 s to a docking-ready folder; Vina (fast preset) 22 s.
