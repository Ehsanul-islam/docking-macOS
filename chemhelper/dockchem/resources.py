"""How many CPU cores the helper may use (default 2, set by the user).

One number governs everything that can run in parallel: Vina's --cpu, OpenMM minimisation threads, RDKit conformer
generation threads, and the thread environment of helper subprocesses (PDB2PQR, Meeko, numpy/BLAS). A request may carry
`cores`; otherwise DOCKCHEM_CORES, otherwise DEFAULT_CORES. Values above the machine's cores are clamped (with a warning);
below 1 or not a whole number is an error.
"""
from __future__ import annotations

import contextvars
import os
import subprocess
import sys

from .errors import InvalidInput

DEFAULT_CORES = 2
_current: contextvars.ContextVar = contextvars.ContextVar("dockchem_cores", default=None)


def total_cores() -> int:
    return os.cpu_count() or 1


def performance_cores() -> int:
    """Performance cores on Apple silicon (efficiency cores slow a docking run down), else all cores."""
    if sys.platform == "darwin":
        try:
            out = subprocess.run(["sysctl", "-n", "hw.perflevel0.physicalcpu"], capture_output=True, text=True, timeout=3)
            if out.returncode == 0 and out.stdout.strip().isdigit():
                return int(out.stdout.strip())
        except Exception:  # noqa: BLE001
            pass
    return total_cores()


def validate(value) -> tuple[int, str | None]:
    """(cores, warning). Raises InvalidInput for anything that is not a whole number >= 1."""
    if isinstance(value, bool) or not isinstance(value, (int, float, str)):
        raise InvalidInput("cores must be a whole number of at least 1")
    try:
        n = float(value)
    except ValueError:
        raise InvalidInput(f"cores must be a whole number of at least 1, not {value!r}") from None
    if n != n or n in (float("inf"), float("-inf")) or n != int(n) or n < 1:
        raise InvalidInput(f"cores must be a whole number of at least 1, not {value!r}")
    n, total = int(n), total_cores()
    if n > total:
        return total, f"cores={n} requested but this computer has {total}; using {total}"
    return n, None


def cores() -> int:
    """Cores for the running job: the request's value, else DOCKCHEM_CORES, else the default of 2."""
    v = _current.get()
    if v is not None:
        return v
    env = os.environ.get("DOCKCHEM_CORES")
    if env:
        try:
            return validate(env)[0]
        except InvalidInput:
            pass
    return min(DEFAULT_CORES, total_cores())


def bind(n: int):
    return _current.set(n)


def unbind(handle) -> None:
    _current.reset(handle)


def thread_env() -> dict:
    n = str(cores())
    return {k: n for k in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS",
                           "NUMEXPR_NUM_THREADS")}


def describe() -> dict:
    total, perf = total_cores(), performance_cores()
    return {"total_cores": total, "performance_cores": perf, "default": min(DEFAULT_CORES, total), "current": cores(),
            "min": 1, "max": total, "advice": "Leave at least one core free so the computer stays responsive; "
                                              "performance cores matter most on Apple silicon."}
