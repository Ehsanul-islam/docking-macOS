#!/usr/bin/env python3
"""Build chemhelper/dockchem/data/ligand_dictionary.json.

An entry is written only if its SMILES matches (heavy-atom composition AND connectivity) a real crystal
ligand in the given PDB files, and the bond-length audit of the result is clean. Each entry records which
structures it was validated on. Candidates that fail are listed and NOT shipped: an unvalidated template
from memory is exactly the kind of silent error this exists to prevent.

  python tools/ligand_dictionary/build_dictionary.py --pdb-dir pdb_structures_50
"""
import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "chemhelper"))
from dockchem import ligand_templates as lt          # noqa: E402
from dockchem.errors import HelperError             # noqa: E402

# resname, SMILES (flat, neutral reference form; stereo comes from the crystal), note, structures to validate on
CANDIDATES = [
    ("STI", "Cc1ccc(NC(=O)c2ccc(CN3CCN(C)CC3)cc2)cc1Nc1nccc(n1)-c1cccnc1", "imatinib, neutral", ["1IEP"]),
    ("AQ4", "COCCOc1cc2ncnc(Nc3cccc(c3)C#C)c2cc1OCCOC", "erlotinib, neutral", ["1M17"]),
    ("BTN", "OC(=O)CCCCC1SCC2NC(=O)NC12", "biotin, neutral acid", ["1STP"]),
    ("BEN", "NC(=N)c1ccccc1", "benzamidine, neutral", ["3PTB"]),
    ("CAU", "CC(C)NCC(O)COc1cccc2[nH]c3ccccc3c12", "carazolol, neutral", ["2RH1"]),
    ("ATP", "Nc1ncnc2c1ncn2C1OC(COP(=O)(O)OP(=O)(O)OP(=O)(O)O)C(O)C1O", "fully protonated phosphates", ["1ATP"]),
    ("GTP", "NC1=NC2=C(N=CN2C2OC(COP(O)(=O)OP(O)(=O)OP(O)(O)=O)C(O)C2O)C(=O)N1", "fully protonated phosphates", ["1JFF", "1TUB"]),
    ("GDP", "NC1=NC2=C(N=CN2C2OC(COP(O)(=O)OP(O)(O)=O)C(O)C2O)C(=O)N1", "fully protonated phosphates", ["1JFF", "1TUB"]),
    ("NAD", "NC(=O)C1=C[N+](=CC=C1)C1OC(COP([O-])(=O)OP(O)(=O)OCC2OC(C(O)C2O)N2C=NC3=C2N=CN=C3N)C(O)C1O", "NAD+, one phosphate anionic", ["4MDH"]),
    ("NDP", "NC(=O)C1=CN(C=CC1)C1OC(COP(O)(=O)OP(O)(=O)OCC2OC(C(OP(O)(O)=O)C2O)N2C=NC3=C2N=CN=C3N)C(O)C1O", "NADPH, protonated phosphates", ["8CAT"]),
    ("AP5", "Nc1ncnc2c1ncn2C1OC(COP(=O)(O)OP(=O)(O)OP(=O)(O)OP(=O)(O)OP(=O)(O)OCC2OC(C(O)C2O)n2cnc3c(N)ncnc23)C(O)C1O", "fully protonated phosphates", ["1AKE"]),
    ("CLR", "CC(C)CCCC(C)C1CCC2C3CC=C4CC(O)CCC4(C)C3CCC12C", "cholesterol", ["2RH1"]),
    ("BU1", "OCCCCO", "1,4-butanediol", ["2RH1"]),
]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pdb-dir", required=True)
    ap.add_argument("--out", default=str(lt.DICT_PATH))
    a = ap.parse_args()
    entries, rejected = {}, []
    for res, smi, note, ids in CANDIDATES:
        ok_on, why = [], []
        for pid in ids:
            f = Path(a.pdb_dir) / f"{pid}.pdb"
            try:
                r = lt.ligand_from_structure(f, res, smiles=smi)
            except HelperError as e:
                why.append(f"{pid}: {e.message[:140]}")
                continue
            if r.audit:
                why.append(f"{pid}: bond-length audit flagged {r.audit[:3]}")
                continue
            ok_on.append(pid)
        if ok_on:
            entries[res] = {"smiles": smi, "protonation": note, "validated_on": ok_on}
        else:
            rejected.append((res, why))
    doc = {"_comment": "Reference templates; each validated against a real crystal ligand (see validated_on). "
                       "Stereochemistry is taken from the structure being prepared, not from here.",
           "version": "2026.1", "entries": dict(sorted(entries.items()))}
    Path(a.out).write_text(json.dumps(doc, indent=2) + "\n")
    print(f"wrote {len(entries)} entries to {a.out}")
    for res, why in rejected:
        print(f"REJECTED {res}: {' | '.join(why)}")


if __name__ == "__main__":
    main()
