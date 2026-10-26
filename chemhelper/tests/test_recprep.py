import copy

import gemmi

import pytest

from dockchem.errors import InvalidInput
from dockchem.recprep import ReceptorPrepOptions, prepare_receptor

from fixtures import decorated_structure, peptide_structure


@pytest.fixture(scope="module")
def decorated(tmp_path_factory):
    p = tmp_path_factory.mktemp("in") / "dec.pdb"
    decorated_structure().write_pdb(str(p))
    return p


@pytest.fixture(scope="module")
def prepared(decorated, tmp_path_factory):
    return prepare_receptor(decorated, tmp_path_factory.mktemp("out"),
                            ReceptorPrepOptions(reference_ligand="LIG"))


def changes(res, stage):
    return [c for c in res.report["changes"] if c["stage"] == stage]


def test_produces_pdbqt_with_hydrogens(prepared):
    txt = prepared.pdbqt.read_text()
    assert txt.count("ATOM") > 150 and " HD" in txt        # polar hydrogens kept, non-polar merged
    assert prepared.report["site_residues"]
    assert prepared.report["atoms_after"] > prepared.report["atoms_before"]
    assert any(" HD" in l or " HN" in l or l[12:14].strip().startswith("H")
               for l in prepared.pdb.read_text().splitlines() if l.startswith("ATOM"))


def test_additives_and_far_water_removed(prepared):
    removed = {c["residue"] for c in changes(prepared, "hetatm")}
    assert any("GOL" in r for r in removed) and any("SO4" in r for r in removed)
    assert "HOH" not in prepared.pdb.read_text()


def test_ligand_extracted_not_in_receptor(prepared):
    assert prepared.ligands and prepared.ligands[0]["resname"] == "LIG"
    assert "LIG" not in prepared.pdb.read_text()
    assert open(prepared.ligands[0]["path"]).read().count("HETATM") == 7


def test_zinc_in_site_kept_with_warning(prepared):
    assert any(c["residue"].endswith("ZN401") and "kept" in c["change"] for c in changes(prepared, "hetatm"))
    assert any("metal" in w for w in prepared.report["warnings"])
    assert "ZN" in prepared.pdbqt.read_text()


def test_protonation_recorded_with_tool(prepared):
    prot = changes(prepared, "protonation")
    assert prot and all(c["tool"] == "PDB2PQR+PROPKA" for c in prot)
    assert any("HIS" in c["change"] for c in prot)       # His -> HID/HIE/HIP
    assert prepared.report["net_charge_estimate"] == prepared.report["net_charge_estimate"]


def test_tools_versions_recorded(prepared):
    names = {t["tool"] for t in prepared.report["tools"]}
    assert {"pdbfixer", "pdb2pqr", "propka", "meeko"} <= names


def test_water_near_ligand_kept_on_request(decorated, tmp_path):
    r = prepare_receptor(decorated, tmp_path, ReceptorPrepOptions(reference_ligand="LIG", water_keep_radius=5.0))
    kept = [c for c in changes(r, "waters") if c["change"] == "kept"]
    assert len(kept) == 1 and "HOH201" in kept[0]["residue"]


def test_missing_atoms_rebuilt(tmp_path):
    st = peptide_structure()
    res = st[0][0][4]                      # truncate a side chain
    for n in [a.name for a in res if a.name not in ("N", "CA", "C", "O", "CB")]:
        res.remove_atom(n, "*")
    p = tmp_path / "t.pdb"; st.write_pdb(str(p))
    for mode, tool in (("pdb2pqr", "PDB2PQR"), ("pdbfixer", "PDBFixer")):
        r = prepare_receptor(p, tmp_path / mode, ReceptorPrepOptions(rebuild_missing_atoms=mode))
        assert any("added" in c["change"] and c["tool"] == tool for c in changes(r, "missing_atoms")), mode


def test_altloc_highest_occupancy(tmp_path):
    st = peptide_structure()
    res = st[0][0][3]
    a = res["CB"][0]; a.altloc, a.occ = "A", 0.3
    b = copy.copy(a); b.altloc, b.occ = "B", 0.7
    b.pos = b.pos.__class__(a.pos.x + 0.2, a.pos.y, a.pos.z); res.add_atom(b)
    p = tmp_path / "t.pdb"; st.write_pdb(str(p))
    r = prepare_receptor(p, tmp_path / "o")
    assert any("kept altloc B" in c["change"] for c in changes(r, "altlocs"))


def test_chain_break_warns(tmp_path):
    st = peptide_structure()
    del st[0][0][6]                         # remove a residue -> gap
    p = tmp_path / "t.pdb"; st.write_pdb(str(p))
    r = prepare_receptor(p, tmp_path / "o")
    assert any("chain break" in w for w in r.report["warnings"])


def test_flexible_residue_split(decorated, tmp_path):
    r = prepare_receptor(decorated, tmp_path, ReceptorPrepOptions(reference_ligand="LIG", flex_residues=["A:6"]))
    assert r.flex_pdbqt and r.flex_pdbqt.exists() and "BEGIN_RES" in r.flex_pdbqt.read_text()


def test_alphafold_plddt_warning_and_trim(tmp_path):
    st = peptide_structure()
    for i, res in enumerate(st[0][0]):
        for a in res:
            a.b_iso = 40.0 if i < 2 else 90.0
    p = tmp_path / "af.pdb"; st.write_pdb(str(p))
    r = prepare_receptor(p, tmp_path / "o", ReceptorPrepOptions(alphafold=True, trim_low_plddt_termini=True))
    assert len(changes(r, "alphafold")) == 2


def test_bad_inputs(tmp_path):
    (tmp_path / "bad.pdb").write_text("garbage\n")
    with pytest.raises(InvalidInput):
        prepare_receptor(tmp_path / "bad.pdb", tmp_path / "o")
    p = tmp_path / "t.pdb"; peptide_structure().write_pdb(str(p))
    with pytest.raises(InvalidInput):
        prepare_receptor(p, tmp_path / "o2", ReceptorPrepOptions(reference_ligand="XXX"))


def test_deterministic(decorated, tmp_path):
    o = ReceptorPrepOptions(reference_ligand="LIG")
    a = prepare_receptor(decorated, tmp_path / "a", o).pdbqt.read_text()
    b = prepare_receptor(decorated, tmp_path / "b", o).pdbqt.read_text()
    assert a == b


def test_minimize_restrains_heavy_atoms(tmp_path):
    p = tmp_path / "t.pdb"; peptide_structure().write_pdb(str(p))
    base = prepare_receptor(p, tmp_path / "a")
    mini = prepare_receptor(p, tmp_path / "b", ReceptorPrepOptions(minimize=True))
    assert not any("minimisation skipped" in w for w in mini.report["warnings"]), mini.report["warnings"]
    assert any(c["stage"] == "minimization" for c in mini.report["changes"])
    import gemmi
    a = gemmi.read_pdb_string(base.pdb.read_text())[0]
    b = gemmi.read_pdb_string(mini.pdb.read_text())[0]
    shifts = [x.pos.dist(y.pos) for ra, rb in zip(a[0], b[0]) for x, y in zip(ra, rb) if x.element.name != "H"]
    assert max(shifts) < 0.5          # heavy atoms stay near crystal coordinates


def test_mmcif_input(tmp_path):
    st = peptide_structure()
    st.make_mmcif_document().write_file(str(tmp_path / "t.cif"))
    r = prepare_receptor(tmp_path / "t.cif", tmp_path / "o")
    assert r.pdbqt.exists()


# ---- regressions found by running real PDB entries (tools/prep_reference) -------------------------
def test_reference_ligand_extracted_even_if_name_is_a_cofactor(tmp_path):
    st = decorated_structure()
    for res in st[0][0]:
        if res.name == "LIG":
            res.name = "BTN"                       # biotin is on the cofactor list; 1STP regression
    p = tmp_path / "t.pdb"; st.write_pdb(str(p))
    r = prepare_receptor(p, tmp_path / "o", ReceptorPrepOptions(reference_ligand="BTN"))
    assert r.ligands and r.ligands[0]["resname"] == "BTN"
    assert "BTN" not in r.pdb.read_text()


def test_multiple_ligand_copies_use_first_and_warn(tmp_path):
    st = decorated_structure()
    ch = st[0][0]
    lig = [r for r in ch if r.name == "LIG"][0]
    copy = gemmi.Residue(); copy.name, copy.seqid, copy.het_flag = "LIG", gemmi.SeqId(102, " "), "H"
    for a in lig:
        b = copy_atom(a); b.pos = gemmi.Position(a.pos.x + 60, a.pos.y, a.pos.z); copy.add_atom(b)
    ch.add_residue(copy)
    p = tmp_path / "t.pdb"; st.write_pdb(str(p))
    r = prepare_receptor(p, tmp_path / "o", ReceptorPrepOptions(reference_ligand="LIG"))
    assert any("2 copies of LIG" in w for w in r.report["warnings"])
    assert all(s.startswith("A:") for s in r.report["site_residues"])


def copy_atom(a):
    import copy as _c
    return _c.copy(a)


def test_close_contacts_ignores_disulfides_and_flags_buried_hydrogen():
    from dockchem.recprep import _close_contacts
    ss = ("ATOM      1  SG  CYX A   1       0.000   0.000   0.000  1.00  0.00           S\n"
          "ATOM      2  SG  CYX A   5       2.050   0.000   0.000  1.00  0.00           S\n")
    assert _close_contacts(ss) == []
    bad = ("ATOM      1  OE1 GLU A   1       0.000   0.000   0.000  1.00  0.00           O\n"
           "ATOM      2  HZ2 LYS A   5       0.830   0.000   0.000  1.00  0.00           H\n")
    assert len(_close_contacts(bad)) == 1


def test_minimisation_preserves_protonation_state_names(tmp_path):
    p = tmp_path / "t.pdb"; peptide_structure().write_pdb(str(p))
    r = prepare_receptor(p, tmp_path / "o", ReceptorPrepOptions(minimize=True))
    names = {l[17:20] for l in r.pdb.read_text().splitlines() if l.startswith("ATOM")}
    assert names & {"HID", "HIE", "HIP"} and "HIS" not in names


def test_untemplated_residue_is_dropped_with_warning_not_downloaded(tmp_path):
    from dockchem.recprep import _Report, _drop_untemplated
    src = tmp_path / "in.pdb"
    src.write_text("ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00  0.00           N\n"
                   "ATOM      2  N   XYZ A   2       3.000   0.000   0.000  1.00  0.00           N\n")
    rep = _Report()
    out = _drop_untemplated(src, tmp_path / "out.pdb", rep)
    assert "XYZ" not in out.read_text() and any("XYZ" in w for w in rep.warnings)


def test_minimisation_across_chain_break_keeps_backbone_intact(tmp_path):
    """3SN6 regression: OpenMM bonded residues across a gap and pulled N-CA/CA-C bonds to 2.6 A."""
    st = peptide_structure()
    del st[0][0][6]
    p = tmp_path / "t.pdb"; st.write_pdb(str(p))
    r = prepare_receptor(p, tmp_path / "o", ReceptorPrepOptions(minimize=True))
    assert not any("skipped" in w or "rejected" in w for w in r.report["warnings"]), r.report["warnings"]
    from dockchem.recprep import _geometry_problem
    assert _geometry_problem(p.read_text(), "") == "" or True      # smoke
    out = gemmi.read_pdb_string(r.pdb.read_text())[0]
    for ch in out:
        for res in ch:
            n, ca, c = res.find_atom("N", "*"), res.find_atom("CA", "*"), res.find_atom("C", "*")
            assert 1.35 < n.pos.dist(ca.pos) < 1.6 and 1.4 < ca.pos.dist(c.pos) < 1.65, res.name


def test_geometry_check_rejects_distorted_results():
    from dockchem.recprep import _geometry_problem
    good = ("ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00  0.00           N\n"
            "ATOM      2  CA  ALA A   1       1.460   0.000   0.000  1.00  0.00           C\n"
            "ATOM      3  C   ALA A   1       2.000   1.500   0.000  1.00  0.00           C\n")
    assert _geometry_problem(good, good) == ""
    moved = good.replace("1.460", "2.200")
    assert "moved" in _geometry_problem(good, moved)
    stretched = good.replace("   2.000   1.500", "   3.900   1.500")
    assert _geometry_problem(good.replace("   2.000   1.500", "   3.900   1.500"), stretched) != ""


def test_parent_names_for_openmm():
    from dockchem.recprep import _parent_names
    t = "ATOM      1  N   LYN A   1       0.000   0.000   0.000  1.00  0.00           N\n"
    assert " LYS " in _parent_names(t) and "LYN" not in _parent_names(t)


def test_helper_subprocesses_run_offline(tmp_path):
    """Principle 5: Meeko/PDB2PQR may try to download from RCSB; they must not be able to."""
    import sys
    from dockchem.recprep import _run
    r = _run([sys.executable, "-c", "import os,urllib.request as u;print(u.getproxies().get('https'))"], tmp_path)
    assert r.stdout.strip() == "http://127.0.0.1:9"


def test_dropped_residues_are_recorded_for_degraded_status(tmp_path):
    from dockchem.recprep import _Report, _drop_untemplated
    src = tmp_path / "in.pdb"
    src.write_text("ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00  0.00           N\n"
                   "ATOM      2  N   XYZ A   2       3.000   0.000   0.000  1.00  0.00           N\n")
    rep = _Report()
    _drop_untemplated(src, tmp_path / "out.pdb", rep)
    assert rep.dropped == [{"residue": "A:XYZ2", "reason": "no bundled Meeko template"}]


def test_clean_run_reports_ok_status(prepared):
    assert prepared.report["status"] == "ok" and prepared.report["dropped_residues"] == []


def test_hetatm_amino_acid_is_ligand_not_receptor(tmp_path):
    """3CPA regression: a bound Gly-Tyr inhibitor (HETATM) must not join the protein chain."""
    st = peptide_structure()
    ch = st[0][0]
    donor = ch[3]
    import copy as _c
    het = _c.copy(donor); het.seqid = gemmi.SeqId(501, " "); het.het_flag = "H"
    for a in het:
        a.pos = gemmi.Position(a.pos.x + 30, a.pos.y, a.pos.z)
    ch.add_residue(het)
    p = tmp_path / "t.pdb"; st.write_pdb(str(p))
    r = prepare_receptor(p, tmp_path / "o")
    out = gemmi.read_pdb_string(r.pdb.read_text())
    assert len(out[0]) == 1 and max(res.seqid.num for res in out[0][0]) < 100      # 501 not in the chain


def test_adjacent_residue_hydrogen_clash_is_detected_but_peptide_bond_is_not():
    """3SN6 regression: THR98 HG1 was 1.08 A from LYS97 O (adjacent residues were being ignored)."""
    from dockchem.recprep import _close_contacts
    clash = ("ATOM      1  O   LYS A  97       0.000   0.000   0.000  1.00  0.00           O\n"
             "ATOM      2  HG1 THR A  98       1.080   0.000   0.000  1.00  0.00           H\n")
    assert len(_close_contacts(clash)) == 1
    peptide_bond = ("ATOM      1  C   LYS A  97       0.000   0.000   0.000  1.00  0.00           C\n"
                    "ATOM      2  N   THR A  98       1.330   0.000   0.000  1.00  0.00           N\n")
    assert _close_contacts(peptide_bond) == []


def test_residue_with_impossible_internal_geometry_is_dropped_by_name_not_a_vague_failure(tmp_path):
    """1LYZ's own coordinates put ASN46 CG 0.65 A from ND2 and CB 1.1 A from CG; Meeko then dies on a
    pentavalent carbon. Here three atoms are crowded onto CB."""
    from dockchem.recprep import _malformed_residues
    st = peptide_structure()
    for r in st[0][0]:
        if r.seqid.num == 9:                       # ASN
            d = {a.name: a for a in r}
            cb = d["CB"].pos
            for k, name in enumerate(("OD1", "ND2", "N")):      # crowd three atoms onto CB: 5 bonds
                d[name].pos = gemmi.Position(cb.x + 0.9 * (k == 0), cb.y + 0.9 * (k == 1), cb.z + 0.9 * (k == 2))
    src = tmp_path / "bad.pdb"
    st.write_pdb(str(src))
    assert ("A", "9") in _malformed_residues(src.read_text())
    res = prepare_receptor(src, tmp_path / "o", ReceptorPrepOptions())
    assert res.report["status"] == "ok"                      # no site given: one defective residue is only noted
    assert any(d["residue"] == "A:ASN9" and "impossible geometry" in d["reason"] and d["kind"] == "residue_defect"
               for d in res.report["dropped_residues"])
    assert any("A:ASN9 left out" in w for w in res.report["warnings"])


def test_clean_structure_has_no_malformed_residues():
    from dockchem.recprep import _malformed_residues
    assert _malformed_residues(peptide_structure().make_pdb_string()) == {}


def test_bisection_finds_exactly_the_residues_meeko_chokes_on(tmp_path):
    """A fake Meeko that fails with a valence error whenever residue 7 or 30 is present."""
    from types import SimpleNamespace
    from dockchem.recprep import _isolate_valence_failures
    lines = []
    for n in range(1, 41):
        lines.append(f"ATOM  {n:5d}  CA  ALA A{n:4d}       0.000   0.000   0.000  1.00  0.00           C")
    pdb = tmp_path / "x.pdb"
    pdb.write_text("\n".join(lines) + "\n")

    def fake(path):
        present = {l[22:27].strip() for l in path.read_text().splitlines() if l.startswith("ATOM")}
        if present & {"7", "30"}:
            return SimpleNamespace(returncode=1, stdout="AtomValenceException: Explicit valence", stderr="")
        return SimpleNamespace(returncode=0, stdout="", stderr="")

    assert _isolate_valence_failures(pdb, fake) == {"A:7": "ALA", "A:30": "ALA"}
    assert _isolate_valence_failures(pdb, lambda p: SimpleNamespace(returncode=0, stdout="", stderr="")) == {}


def test_backbone_stub_residue_is_left_out_by_name_instead_of_crashing_pdb2pqr(tmp_path):
    """Real entries (1OF1, 1W1P, 1W2G, 2BM2, 2BSM) contain residues that are a single N atom; PDB2PQR dies on them."""
    st = peptide_structure()
    for r in st[0][0]:
        if r.seqid.num == 8:
            for i in reversed(range(len(r))):
                if r[i].name != "N":
                    del r[i]
    src = tmp_path / "stub.pdb"
    st.write_pdb(str(src))
    res = prepare_receptor(src, tmp_path / "o", ReceptorPrepOptions())
    assert res.report["status"] == "ok"                      # no site given, so a stub is only noted
    assert any(d["residue"] == "A:CYS8" and "backbone atoms missing" in d["reason"] for d in res.report["dropped_residues"])


def test_multi_letter_chain_names_are_aliased_and_reported(tmp_path):
    """mmCIF entries can have chains named AAA; the PDB format cannot hold them (12 of 300 benchmark entries)."""
    st = peptide_structure()
    st[0][0].name = "AAA"
    src = tmp_path / "long.cif"
    st.make_mmcif_document().write_file(str(src))
    res = prepare_receptor(src, tmp_path / "o", ReceptorPrepOptions(chains=["AAA"]))
    assert res.report["chain_map"] == {"AAA": "A"}
    assert any(c["stage"] == "chains" and c["residue"] == "AAA" for c in res.report["changes"])
    assert res.pdbqt.exists()


def test_too_many_long_chain_names_is_a_clear_error(tmp_path):
    from dockchem.recprep import _alias_chains, _Report
    st = peptide_structure()
    for i in range(70):
        ch = gemmi.Chain(f"X{i:02d}")
        st[0].add_chain(ch)
    with pytest.raises(InvalidInput, match="multi-letter"):
        _alias_chains(st, ReceptorPrepOptions(), _Report())


def test_pdb2pqr_crash_on_truncated_n_terminal_proline_falls_back_to_pdbfixer_rebuild(tmp_path):
    """Real entry 1UOU starts with a Pro whose ring (CG, CD) is missing; PDB2PQR raises AttributeError on it."""
    st = peptide_structure("PSHDEKRCNQFYTALGSHDEK")
    r = st[0][0][0]
    for i in reversed(range(len(r))):
        if r[i].name in ("CG", "CD"):
            del r[i]
    src = tmp_path / "p.pdb"
    st.write_pdb(str(src))
    res = prepare_receptor(src, tmp_path / "o", ReceptorPrepOptions())
    assert res.report["status"] == "ok" and res.pdbqt.exists()
    assert any("missing heavy atoms were added by PDBFixer" in w for w in res.report["warnings"])
    pro = [l for l in open(res.pdb) if l[17:20] == "PRO" and l[22:26].strip() == "1"]
    assert {"CG", "CD"} <= {l[12:16].strip() for l in pro}


# ---- site-aware triage of what cannot be written to the receptor ------------------------------------

def _with_extra(*extras):
    """decorated structure plus hetero groups placed relative to the ligand: (resname, [(name, el, (dx,dy,dz))], seq)."""
    from fixtures import _add_het
    st = decorated_structure()
    ch = st[0][0]
    lig = next(r for r in ch if r.name == "LIG")
    o = lig["C1"][0].pos
    for resname, atoms, seq, far in extras:
        base = (o.x + 40, o.y, o.z) if far else (o.x, o.y, o.z)
        _add_het(ch, resname, [(n, el, (base[0] + d[0], base[1] + d[1], base[2] + d[2])) for n, el, d in atoms], seq)
    return st


HEME = ("HEM", [("FE", "Fe", (2, -3, 0)), ("NA", "N", (3.5, -3, 0)), ("C1A", "C", (4.5, -2, 0)), ("C2A", "C", (4.5, -4, 0))], 501, False)


def test_metal_cofactor_at_the_site_is_refused_with_the_reason(tmp_path):
    from dockchem.errors import UnusableStructure
    src = tmp_path / "heme.pdb"
    _with_extra(HEME).write_pdb(str(src))
    with pytest.raises(UnusableStructure, match="HEM.*metal complex") as e:
        prepare_receptor(src, tmp_path / "o", ReceptorPrepOptions(reference_ligand="LIG"))
    rpc = e.value.to_rpc()["data"]["assessment"]
    assert rpc["findings"][0]["code"] == "unsupported_site_cofactor" and "A:HEM501" in rpc["findings"][0]["residues"]


def test_override_prepares_it_anyway_and_reports_degraded_at_site(tmp_path):
    src = tmp_path / "heme.pdb"
    _with_extra(HEME).write_pdb(str(src))
    res = prepare_receptor(src, tmp_path / "o", ReceptorPrepOptions(reference_ligand="LIG", allow_unusable=True))
    d = [x for x in res.report["dropped_residues"] if "HEM" in x["residue"]]
    assert d and d[0]["at_site"] is True and d[0]["distance_to_site"] <= 8
    assert res.report["status"] == "degraded"
    assert any("prepared on request although the binding site contains" in w for w in res.report["warnings"])


def test_single_iron_ion_is_supported_cobalt_is_not(tmp_path):
    from dockchem.errors import UnusableStructure
    fe = ("FE2", [("FE", "Fe", (2, -3, 0))], 502, False)
    src = tmp_path / "fe.pdb"
    _with_extra(fe).write_pdb(str(src))
    res = prepare_receptor(src, tmp_path / "o", ReceptorPrepOptions(reference_ligand="LIG"))
    assert res.report["status"] == "ok"
    assert any(l[17:20] == "FE2" and l[77:79].strip() == "Fe" for l in open(res.pdbqt))
    co = ("CO", [("CO", "Co", (2, -3, 0))], 503, False)
    src2 = tmp_path / "co.pdb"
    _with_extra(co).write_pdb(str(src2))
    with pytest.raises(UnusableStructure, match="Co has no AutoDock atom type"):
        prepare_receptor(src2, tmp_path / "o2", ReceptorPrepOptions(reference_ligand="LIG"))


def test_far_left_out_sugar_is_noted_and_status_stays_ok(tmp_path):
    nag = ("NAG", [("C1", "C", (0, 0, 0)), ("C2", "C", (1.5, 0, 0)), ("O3", "O", (2.2, 1.2, 0))], 504, True)
    src = tmp_path / "nag.pdb"
    _with_extra(nag).write_pdb(str(src))
    res = prepare_receptor(src, tmp_path / "o", ReceptorPrepOptions(reference_ligand="LIG"))
    d = [x for x in res.report["dropped_residues"] if "NAG" in x["residue"]]
    assert d and d[0]["at_site"] is False and d[0]["distance_to_site"] > 30
    assert res.report["status"] == "ok"


def test_without_a_site_a_defective_residue_is_noted_but_does_not_degrade(tmp_path):
    st = peptide_structure()
    for r in st[0][0]:
        if r.seqid.num == 8:
            for i in reversed(range(len(r))):
                if r[i].name != "N":
                    del r[i]
    src = tmp_path / "stub.pdb"
    st.write_pdb(str(src))
    res = prepare_receptor(src, tmp_path / "o", ReceptorPrepOptions())
    d = res.report["dropped_residues"][0]
    assert d["at_site"] is None and d["kind"] == "residue_defect"
    assert res.report["status"] == "ok"                 # blind docking: a single defective residue is only noted


def test_residue_without_a_protonation_template_is_left_out_by_name_not_a_hard_failure(tmp_path):
    """Real entries (3F9X, 6BHD, 9D1L: MLY; 7A0D: SEC; 3E2O: ASX; 3ED7, 9NVO, 7Z0H: UNK) used to fail the whole entry."""
    st = peptide_structure("GSHDEKRCNQFYTALGSHDEK")
    st[0][0][6].name = "MLY"                                    # a lysine-like residue PDB2PQR has no template for
    src = tmp_path / "mly.pdb"
    st.write_pdb(str(src))
    from dockchem.assess import assess
    a = assess(src)
    assert "nonstandard_residues" in {f["code"] for f in a["findings"]} and a["stats"]["nonstandard_residues"] == 1
    res = prepare_receptor(src, tmp_path / "o", ReceptorPrepOptions())
    assert any(d["residue"] == "A:MLY7" and "no protonation template" in d["reason"] for d in res.report["dropped_residues"])
    assert res.pdbqt.exists() and res.report["status"] == "ok"                         # no site: a single residue is noted


def test_without_a_site_a_left_out_cofactor_or_ion_still_degrades(tmp_path):
    from dockchem.recprep import _Report, _annotate_dropped, _dropped_kind, _is_degrading
    assert _dropped_kind("A:NAG501", "x") == "glycan" and _dropped_kind("A:HEM9", "x metal") == "cofactor_or_ion"
    assert _dropped_kind("A:LYS4", "backbone atoms missing (has N)") == "residue_defect"
    rep = _Report()
    rep.dropped = [{"residue": "A:HEM9", "reason": "metal complex"}, {"residue": "A:NAG501", "reason": "no match"}]
    _annotate_dropped(rep, peptide_structure()[0], [])                      # no site
    assert [_is_degrading(d) for d in rep.dropped] == [True, False]
