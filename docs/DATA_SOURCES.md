# Bundled data sources

| File | What | Source | Licence | Refresh |
|---|---|---|---|---|
| `chemhelper/dockchem/data/ccd_ligands.json.gz` | 51,493 ligand definitions (flat SMILES, formula, name) used to establish bond orders for ligands taken from PDB/mmCIF files | wwPDB Chemical Component Dictionary, `https://files.wwpdb.org/pub/pdb/data/monomers/components.cif.gz`, file dated 2026-10-17 | wwPDB archive data are released under CC0 1.0 (public domain). **[VERIFY with counsel before commercial release]** | `python tools/ligand_dictionary/build_ccd.py --ccd components.cif.gz --last-modified "<header>"` |

Residue templates for the receptor come from Meeko's bundled `residue_chem_templates.json` (see its licence in
plan §16).
