import json
import os
import stat
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))
import run_equivalence as eq


def fake_vina(path, score):
    path.write_text(f"""#!/usr/bin/env python3
import sys
a = sys.argv
out = a[a.index("--out") + 1]
open(out, "w").write("REMARK VINA RESULT:    {score}      0.000      0.000\\nATOM\\n")
""")
    path.chmod(path.stat().st_mode | stat.S_IEXEC)
    return str(path)


def make_cases(tmp):
    (tmp / "cases.json").write_text(json.dumps({"cases": [{
        "name": "c1", "receptor": "r.pdbqt", "ligand": "l.pdbqt",
        "center": [0, 1, 2], "size": [20, 20, 20], "seed": 1}]}))
    return str(tmp / "cases.json")


def test_exact_pass(tmp_path):
    a, b = fake_vina(tmp_path / "a", "-7.500"), fake_vina(tmp_path / "b", "-7.500")
    assert eq.main(["--reference", a, "--candidate", b, "--cases", make_cases(tmp_path)]) == 0


def test_exact_fails_on_any_difference(tmp_path):
    a, b = fake_vina(tmp_path / "a", "-7.500"), fake_vina(tmp_path / "b", "-7.501")
    assert eq.main(["--reference", a, "--candidate", b, "--cases", make_cases(tmp_path)]) == 1


def test_tolerance_mode(tmp_path):
    a, b = fake_vina(tmp_path / "a", "-7.500"), fake_vina(tmp_path / "b", "-7.560")
    cases = make_cases(tmp_path)
    assert eq.main(["--reference", a, "--candidate", b, "--cases", cases, "--mode", "tolerance"]) == 0
    assert eq.main(["--reference", a, "--candidate", b, "--cases", cases,
                    "--mode", "tolerance", "--tol", "0.01"]) == 1
