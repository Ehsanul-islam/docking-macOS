"""Malformed input must never produce an InternalError (a bug), a crash or a hang: only results or typed errors."""
import random
import time

import pytest

from dockchem.server import Server
from fixtures import peptide_structure

OK_ERRORS = {"InvalidInput", "UnusableStructure", "FormatLimit", "ToolFailure", "SelectionError", "NeedsTemplate",
             "UnsupportedFeature", "Cancelled"}


def call(method, **params):
    t0 = time.time()
    r = Server().handle({"jsonrpc": "2.0", "id": 1, "method": method, "params": params})
    assert time.time() - t0 < 60, f"{method} took {time.time() - t0:.0f}s"
    assert r is not None and ("result" in r or "error" in r)
    if "error" in r:
        data = r["error"].get("data") or {}
        assert data.get("type") in OK_ERRORS or r["error"]["code"] in (-32601, -32700), r["error"]
    return r


@pytest.fixture(scope="module")
def good_pdb():
    return peptide_structure("GSHDEKRCNQ").make_pdb_string()


def corpus(good):
    lines = good.splitlines()
    atoms = [l for l in lines if l.startswith("ATOM")]
    yield "empty", b""
    yield "whitespace", b"  \n\n \n"
    yield "binary", bytes(range(256)) * 20
    yield "nulls", b"\x00" * 500
    yield "one_line", atoms[0].encode()
    yield "no_coords", b"ATOM      1  CA  ALA A   1\nEND\n"
    yield "nan_coords", "\n".join(l[:30] + "     nan     nan     nan" + l[54:] for l in atoms).encode()
    yield "huge_coords", "\n".join(l[:30] + f"{1e7:8.0f}{-1e7:8.0f}{1e7:8.0f}"[:24] + l[54:] for l in atoms).encode()
    yield "all_same_xyz", "\n".join(l[:30] + "   0.000   0.000   0.000" + l[54:] for l in atoms).encode()
    yield "crlf", good.replace("\n", "\r\n").encode()
    yield "unicode", ("REMARK " + "é漢字🧪" * 50 + "\n" + good).encode()
    yield "very_long_line", ("REMARK " + "x" * 200000 + "\n" + good).encode()
    yield "truncated_mid_line", good.encode()[: len(good) // 2 + 7]
    yield "only_hetatm_waters", "\n".join(f"HETATM{i:5d}  O   HOH A{i:4d}    {i:8.3f}   0.000   0.000  1.00  0.00           O" for i in range(1, 30)).encode()
    yield "unknown_residues", "\n".join(l[:17] + "ZZZ" + l[20:] for l in atoms).encode()
    yield "duplicate_atoms", "\n".join(atoms + atoms).encode()
    yield "no_element_column", "\n".join(l[:76] for l in atoms).encode()
    yield "cif_garbage", b"data_x\nloop_\n_atom_site.id\n1 2 3\n"
    yield "html", b"<html><body>404 not found</body></html>"


def test_malformed_structure_files_never_cause_internal_errors(tmp_path, good_pdb):
    for name, data in corpus(good_pdb):
        for ext in (".pdb", ".cif"):
            f = tmp_path / f"{name}{ext}"
            f.write_bytes(data)
            for m in ("structure.assess", "format.sniff"):
                call(m, path=str(f))
            call("format.convert", source=str(f), dest=str(tmp_path / f"{name}{ext}.out.sdf"))
            call("ligand.prepare", path=str(f))
            call("library.prepare", path=str(f))


def test_random_mutations_of_a_valid_structure_never_cause_internal_errors(tmp_path, good_pdb):
    rng = random.Random(7)
    raw = good_pdb.encode()
    for i in range(60):
        b = bytearray(raw)
        for _ in range(rng.randint(1, 8)):
            op = rng.choice(["flip", "delete", "dup", "truncate", "nan", "insert"])
            pos = rng.randrange(len(b))
            if op == "flip":
                b[pos] = rng.randrange(256)
            elif op == "delete":
                del b[pos:pos + rng.randint(1, 120)]
            elif op == "dup":
                b[pos:pos] = b[pos:pos + rng.randint(1, 200)]
            elif op == "truncate":
                b = b[:max(10, pos)]
            elif op == "nan":
                b[pos:pos] = b" nan "
            else:
                b[pos:pos] = bytes(rng.randrange(256) for _ in range(rng.randint(1, 30)))
            if not b:
                b = bytearray(b"\n")
        f = tmp_path / f"mut{i}.pdb"
        f.write_bytes(bytes(b))
        call("structure.assess", path=str(f))
        call("format.sniff", path=str(f))


@pytest.mark.parametrize("smiles", ["", "(((", "C" * 5000, "C\x00C", "[Xx]", "C1CC", "c1ccccc1" * 200, "O=C=O.[Na+].[Cl-]",
                                    "[H][H]", "*", "C%99", "CC(C)(C)(C)(C)C", "\ud800", "[U]", "N#N", "ClCl"])
def test_malformed_smiles_never_cause_internal_errors(smiles):
    call("ligand.inspect", smiles=smiles)
    call("ligand.prepare", smiles=smiles)


@pytest.mark.parametrize("expr", ["", "(", "all and", "resi 1-", "within -3 of all", "chain " + "A," * 3000, "not not not all",
                                  "resn \x00", "bfactor > nan", "((((((((all))))))))", "name CA or " * 400 + "all", "resi 999999999999999999999"])
def test_malformed_selections_are_selection_errors_not_bugs(expr):
    call("selection.validate", expr=expr)


def test_bad_parameter_types_are_typed_errors(tmp_path):
    for params in ({"path": 5}, {"path": ["a"]}, {"path": None}, {"path": {"x": 1}}):
        call("structure.assess", **params)
        call("format.sniff", **params)
    call("ligand.prepare", smiles=12345)
    call("ligand.inspect", smiles=["C"])
    call("receptor.prepare", path=str(tmp_path), outdir=str(tmp_path / "o"), options="nope")
    call("receptor.prepare", path=str(tmp_path), outdir=str(tmp_path / "o"), options={"ph": "high"})


def test_receptor_prepare_survives_a_few_hostile_but_parseable_files(tmp_path, good_pdb):
    for name in ("all_same_xyz", "unknown_residues", "duplicate_atoms", "no_element_column"):
        data = dict(corpus(good_pdb))[name]
        f = tmp_path / f"{name}.pdb"
        f.write_bytes(data)
        call("receptor.prepare", path=str(f), outdir=str(tmp_path / f"o_{name}"))
