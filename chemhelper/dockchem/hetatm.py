"""HETATM classification (plan §8.1).

Curated lists, versioned: changing them is a scientific-default change (plan §12.4) and needs a
benchmark re-run. The plan calls for the full PDB Chemical Component Dictionary; until a CCD
snapshot is bundled (open item) this module falls back to heavy-atom count for 'ligand'.
"""
LIST_VERSION = "2026.1"

WATER = {"HOH", "WAT", "DOD", "H2O", "TIP", "TIP3"}
ADDITIVES = {"SO4", "PO4", "GOL", "EDO", "PEG", "PGE", "PG4", "ACT", "DMS", "MPD", "CL", "BME",
             "TRS", "FMT", "ACY", "IPA", "EOH", "NO3", "SCN", "IOD", "BR", "NA", "K", "PGO", "1PE",
             "P6G", "MES", "EPE", "CIT", "TAR", "MLI", "BU3", "OXL"}
METALS = {"ZN", "MG", "MN", "FE", "FE2", "CA", "CU", "CU1", "CO", "NI", "CD", "HG", "MO", "W", "PT"}
COFACTORS = {"HEM", "HEC", "NAD", "NAP", "NDP", "NAI", "FAD", "FMN", "ATP", "ADP", "AMP", "ANP",
             "GTP", "GDP", "SAM", "SAH", "COA", "PLP", "TPP", "BTN", "FES", "SF4", "F3S", "CLA",
             "BCL", "THF", "H4B"}
GLYCANS = {"NAG", "MAN", "BMA", "FUC", "GAL", "GLC", "SIA", "NDG", "BGC", "XYP", "A2G"}
# Modified residues kept as polymer after mapping to the parent (PDBFixer performs the mapping).
MODIFIED = {"MSE": "MET", "SEP": "SER", "TPO": "THR", "PTR": "TYR", "CSO": "CYS", "HYP": "PRO",
            "MLY": "LYS", "KCX": "LYS", "CME": "CYS", "SMC": "CYS", "FME": "MET", "ALY": "LYS"}

MIN_LIGAND_HEAVY = 6


def classify(resname: str, n_heavy: int, overrides: dict | None = None) -> str:
    """Return one of: water, additive, metal, cofactor, glycan, modified, ligand, small_unknown.
    `overrides` maps residue name -> class (users can override any classification)."""
    r = resname.upper().strip()
    if overrides and r in overrides:
        return overrides[r]
    if r in WATER:
        return "water"
    if r in MODIFIED:
        return "modified"
    if r in METALS:
        return "metal"
    if r in COFACTORS:
        return "cofactor"
    if r in GLYCANS:
        return "glycan"
    if r in ADDITIVES:
        return "additive"
    return "ligand" if n_heavy >= MIN_LIGAND_HEAVY else "small_unknown"
