"""Provenance records (plan §14)."""
import hashlib
import platform
import time
from pathlib import Path

from . import __version__


def sha256_file(path: str | Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def record(tool: str, version: str, params: dict | None = None, **extra) -> dict:
    return {
        "tool": tool,
        "version": version,
        "params": params or {},
        "app_version": __version__,
        "arch": platform.machine(),
        "timestamp": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        **extra,
    }
