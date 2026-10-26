"""Cooperative cancellation. The server gives every running job a token; long operations call check() at safe
points (stage boundaries, minimisation chunks) and poll it while waiting for helper subprocesses."""
from __future__ import annotations

import contextvars
import threading

from .errors import Cancelled


class Token:
    def __init__(self):
        self._event = threading.Event()

    def cancel(self) -> None:
        self._event.set()

    @property
    def cancelled(self) -> bool:
        return self._event.is_set()

    def raise_if_cancelled(self) -> None:
        if self._event.is_set():
            raise Cancelled("job cancelled")


_current: contextvars.ContextVar = contextvars.ContextVar("dockchem_cancel_token", default=None)


def current() -> Token | None:
    return _current.get()


def check() -> None:
    """Raise Cancelled if the running job was cancelled. A no-op outside a server job."""
    t = _current.get()
    if t is not None:
        t.raise_if_cancelled()


def bind(token: Token | None):
    """Make `token` the current job's token; returns a reset handle for unbind()."""
    return _current.set(token)


def unbind(handle) -> None:
    _current.reset(handle)
