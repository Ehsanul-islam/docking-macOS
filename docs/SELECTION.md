# Selection and editing

One general way to select and delete parts of a structure (chains, hetero atoms, residue names or ranges,
elements, "within X Å of a ligand"), used by the app's selection UI and by receptor preparation.

## Selection language

Keywords are case-insensitive. `and` binds tighter than `or`; use parentheses to be explicit.

| Form | Meaning |
|---|---|
| `all`, `none` | every atom / no atom |
| `chain A,B` | chains |
| `resn HOH,GOL,SO*` | residue names (`*` and `?` wildcards) |
| `resi 10-20,45,52A` | residue numbers, ranges, insertion codes |
| `name CA,CB`, `elem C,N`, `altloc A` | atom name, element, alternate location |
| `protein` `nucleic` `water` `hetero` | broad classes (`hetero` = HETATM records) |
| `ligand` `cofactor` `metal` `additive` `glycan` `modified` | classes from `hetatm.py` (versioned lists) |
| `backbone` `sidechain` `hydrogen` | protein backbone / side-chain heavy atoms / hydrogens |
| `bfactor > 60`, `occ < 1` | numeric tests (`<` `>` `<=` `>=`) |
| `within 5 of <sel>` | atoms within 5 Å of the selection (includes the selection) |
| `around 5 of <sel>` | the same, excluding the selection |
| `byres <sel>` | whole residues containing any selected atom |
| `not <sel>`, `<a> and <b>`, `<a> or <b>`, `( <sel> )` | boolean logic |

Examples: `delete chain B` · `delete hetero and not within 5 of resn STI` ·
`delete water and not within 4 of ligand` · `keep byres within 8 of resn BTN` · `delete resi 1-10`.

## Edits, preview, undo

An edit is `delete <sel>` or `keep <sel>`. Editing is **non-destructive**: the list of edit commands is the entire
state. Undo drops the last command, redo puts it back, and replaying the list always yields the same structure.
The input file is never modified.

`selection.preview` returns the counts, the affected residue ids (for highlighting in the 3D view) and warnings:
partial residues, new chain breaks, removal of pocket residues (within 5 Å of a ligand), removal of most or all of
the protein. `selection.apply` writes the edited structure. Both are stateless: the app sends its edit list.

In receptor preparation the same list is the `edits` option; each edit is logged in the Prep Report
(`stage: selection`) with its warnings.

## Limits

- Selection is by atom/residue properties and distance; there is no sequence-pattern or secondary-structure term.
- Distances use the coordinates in the file (first altloc already chosen upstream); symmetry mates are not considered.
- `partial` residues (some atoms removed) are warned about but allowed, since a user may want them; receptor
  preparation will usually fail on such residues.
- The 3D highlighting and the sequence strip are viewer features and are not built yet.
