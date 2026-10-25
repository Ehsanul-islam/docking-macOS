#!/usr/bin/env python3
"""Build a diverse protein test set from the PDB itself: the most-studied proteins per organism group.

"Most used" here means: proteins with the most PDB structures (a UniProt group's size), one representative per
protein (best resolution). Per organism group (human, other animals, bacteria, plants, others) and per method
(X-ray, cryo-EM, NMR) the largest groups are taken first.

  python tools/batch/make_manifest_rcsb.py --out tools/batch/sets/most_studied_1000.json --n 1000
"""
import argparse
import json
import urllib.request

API = "https://search.rcsb.org/rcsbsearch/v2/query"
LINEAGE = "rcsb_entity_source_organism.taxonomy_lineage.id"

# group -> (quota share, taxonomy include ids, taxonomy exclude ids)
GROUPS = {
    "human": (0.40, ["9606"], []),
    "animal": (0.15, ["33208"], ["9606"]),          # Metazoa excluding human
    "bacteria": (0.25, ["2"], []),
    "plant": (0.05, ["33090"], []),                 # Viridiplantae
    "other": (0.15, ["10239", "4751", "2157"], []),  # viruses, fungi, archaea
}
METHODS = {                                          # share, exptl.method, max resolution, needs ligand
    "xray": (0.90, "X-RAY DIFFRACTION", 2.5, True),
    "cryoem": (0.07, "ELECTRON MICROSCOPY", 3.5, True),
    "nmr": (0.03, "SOLUTION NMR", None, False),
}
MAX_ATOMS = 100_000                                  # bounds run time; recorded as a limit of this run


def term(attr, op, value, negate=False):
    p = {"attribute": attr, "operator": op, "value": value}
    if negate:
        p["negation"] = True
    return {"type": "terminal", "service": "text", "parameters": p}


def query(inc, exc, method, maxres, ligand, start, rows):
    nodes = [term("exptl.method", "exact_match", method),
             term("rcsb_entry_info.deposited_atom_count", "less_or_equal", MAX_ATOMS)]
    nodes.append({"type": "group", "logical_operator": "or", "nodes": [term(LINEAGE, "exact_match", t) for t in inc]}
                 if len(inc) > 1 else term(LINEAGE, "exact_match", inc[0]))
    for t in exc:
        nodes.append(term(LINEAGE, "exact_match", t, negate=True))
    if maxres:
        nodes.append(term("rcsb_entry_info.resolution_combined", "less_or_equal", maxres))
    if ligand:
        nodes.append(term("rcsb_entry_info.nonpolymer_entity_count", "greater_or_equal", 1))
    return {"query": {"type": "group", "logical_operator": "and", "nodes": nodes}, "return_type": "polymer_entity",
            "request_options": {"group_by": {"aggregation_method": "matching_uniprot_accession",
                                             "ranking_criteria_type": {"sort_by": "rcsb_entry_info.resolution_combined",
                                                                        "direction": "asc"}},
                                "group_by_return_type": "groups", "paginate": {"start": start, "rows": rows}}}


def post(q):
    req = urllib.request.Request(API, data=json.dumps(q).encode(), headers={"Content-Type": "application/json"})
    with urllib.request.urlopen(req, timeout=120) as r:
        return json.load(r)


def groups(inc, exc, method, maxres, ligand):
    out, start = [], 0
    while True:
        r = post(query(inc, exc, method, maxres, ligand, start, 1000))
        gs = r.get("group_set", [])
        out += [(g["count"], g["identifier"], g["representative"]["identifier"]) for g in gs]
        start += 1000
        if not gs or start >= r.get("group_by_count", 0):
            return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("--n", type=int, default=1000)
    a = ap.parse_args()
    items, seen = [], set()
    for gname, (gshare, inc, exc) in GROUPS.items():
        for mname, (mshare, method, maxres, ligand) in METHODS.items():
            quota = round(a.n * gshare * mshare)
            if not quota:
                continue
            gs = sorted(groups(inc, exc, method, maxres, ligand), reverse=True)
            taken = 0
            for count, uniprot, rep in gs:
                pdb = rep.split("_")[0]
                if pdb in seen:
                    continue
                seen.add(pdb)
                items.append({"id": pdb, "organism_group": gname, "method": mname, "uniprot": uniprot,
                              "structures_for_protein": count})
                taken += 1
                if taken >= quota:
                    break
            print(f"{gname:9s}{mname:7s} quota {quota:4d} taken {taken:4d} (of {len(gs)} proteins)")
    json.dump({"description": "most-studied proteins per organism group and method (UniProt group size), one best-resolution "
                              f"representative each; entries over {MAX_ATOMS} atoms excluded to bound run time",
               "source": "RCSB PDB search API", "items": items}, open(a.out, "w"), indent=1)
    print(len(items), "entries ->", a.out)


if __name__ == "__main__":
    main()
