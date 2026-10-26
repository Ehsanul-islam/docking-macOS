"""Input guards: reject what cannot be processed with a clear message instead of failing deep inside a tool."""
from __future__ import annotations

import os
from pathlib import Path

from .errors import InvalidInput

MAX_INPUT_BYTES = int(os.environ.get("DOCKCHEM_MAX_INPUT_MB", "300")) * 1024 * 1024
STAGE_TIMEOUT_S = float(os.environ.get("DOCKCHEM_STAGE_TIMEOUT", "1800"))     # per helper subprocess


def check_input_file(path, what: str = "input file") -> Path:
    """Exists, is a regular readable file, is not absurdly large. Returns the Path."""
    if path is None or str(path).strip() == "":
        raise InvalidInput(f"{what}: no path given")
    p = Path(str(path))
    try:
        if not p.exists():
            raise InvalidInput(f"{what} not found: {p}")
        if p.is_dir():
            raise InvalidInput(f"{what} is a folder, not a file: {p}")
        size = p.stat().st_size
        with open(p, "rb") as fh:
            fh.read(1)
    except InvalidInput:
        raise
    except PermissionError as e:
        raise InvalidInput(f"{what} cannot be read (permission denied): {p}") from e
    except OSError as e:
        raise InvalidInput(f"{what} cannot be read: {p} ({e.strerror or e})") from e
    if size == 0:
        raise InvalidInput(f"{what} is empty: {p}")
    if size > MAX_INPUT_BYTES:
        raise InvalidInput(f"{what} is {size / 1e6:.0f} MB, over the {MAX_INPUT_BYTES / 1e6:.0f} MB limit: {p}")
    return p


def check_structure_text(path, what: str = "structure file") -> None:
    """PDB-style coordinate records must be plain text; stray bytes make the parser fail late and obscurely."""
    p = Path(str(path))
    if p.suffix.lower() not in (".pdb", ".ent", ".pqr", ".pdbqt"):
        return
    with open(p, "rb") as fh:
        for n, line in enumerate(fh, 1):
            if line.startswith((b"ATOM", b"HETATM", b"ANISOU")) and any(b > 127 or b == 0 for b in line):
                raise InvalidInput(f"{what} has non-text bytes in a coordinate record (line {n}); it is damaged or not a PDB file")
