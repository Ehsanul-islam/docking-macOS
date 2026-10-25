"""Generic resumable batch runner with hard per-item timeouts.

Each item is processed by `fn(item) -> dict` in a worker process. A worker that exceeds `timeout` seconds,
crashes (segfault, OOM kill) or raises is recorded as a result of its own (category timeout / crash /
exception) and replaced; the batch never stops for one bad item. Results are appended to a JSONL file
(one line per item id), so an interrupted run resumes where it left off.

fn must be importable at module level in the worker (multiprocessing 'spawn' safe).
"""
from __future__ import annotations

import json
import multiprocessing as mp
import os
import time
import traceback
from pathlib import Path


def _worker(fn_path, conn):
    mod, _, name = fn_path.rpartition(":")
    import importlib
    fn = getattr(importlib.import_module(mod), name)
    while True:
        try:
            item = conn.recv()
        except EOFError:
            return
        if item is None:
            return
        t0 = time.time()
        try:
            res = fn(item)
            res.setdefault("id", item["id"])
        except Exception as e:  # noqa: BLE001
            res = {"id": item["id"], "category": "exception", "error": f"{type(e).__name__}: {e}",
                   "trace": traceback.format_exc()[-800:]}
        res["seconds"] = round(time.time() - t0, 2)
        conn.send(res)


class _Slot:
    def __init__(self, fn_path):
        self.fn_path = fn_path
        self.start()

    def start(self):
        self.parent, child = mp.Pipe()
        self.proc = mp.get_context("spawn").Process(target=_worker, args=(self.fn_path, child), daemon=True)
        self.proc.start()
        self.item = None
        self.t0 = 0.0

    def assign(self, item):
        self.item, self.t0 = item, time.time()
        self.parent.send(item)

    def kill(self):
        try:
            self.proc.kill()
            self.proc.join(5)
        finally:
            self.parent.close()


def run_batch(items: list[dict], fn_path: str, out: Path, workers: int = 4, timeout: float = 120,
              progress_every: int = 200, log=print) -> dict:
    """items: [{'id': ..., ...}]; fn_path: 'module:function'. Returns {'done', 'skipped', 'seconds'}."""
    out = Path(out)
    done_ids = set()
    if out.exists():
        for line in out.read_text().splitlines():
            try:
                done_ids.add(json.loads(line)["id"])
            except Exception:  # noqa: BLE001  truncated last line from an interrupted run
                pass
    todo = [i for i in items if i["id"] not in done_ids]
    log(f"{len(items)} items, {len(done_ids)} already done, {len(todo)} to run, {workers} workers")
    if not todo:
        return {"done": 0, "skipped": len(done_ids), "seconds": 0}
    t_start, n_done = time.time(), 0
    slots = [_Slot(fn_path) for _ in range(min(workers, len(todo)))]
    queue = list(reversed(todo))
    with out.open("a") as fh:
        def write(res):
            nonlocal n_done
            fh.write(json.dumps(res) + "\n")
            fh.flush()
            n_done += 1
            if n_done % progress_every == 0:
                rate = n_done / (time.time() - t_start)
                log(f"  {n_done}/{len(todo)}  {rate:.1f}/s  eta {(len(todo) - n_done) / rate / 60:.0f} min")
        for s in slots:
            if queue:
                s.assign(queue.pop())
        while any(s.item is not None for s in slots):
            for s in slots:
                if s.item is None:
                    continue
                if s.parent.poll(0.05):
                    try:
                        write(s.parent.recv())
                    except EOFError:
                        write({"id": s.item["id"], "category": "crash", "error": "worker died"})
                        s.kill(); s.start()
                    s.item = None
                elif not s.proc.is_alive():
                    write({"id": s.item["id"], "category": "crash",
                           "error": f"worker exited with code {s.proc.exitcode}"})
                    s.kill(); s.start(); s.item = None
                elif time.time() - s.t0 > timeout:
                    write({"id": s.item["id"], "category": "timeout", "error": f"exceeded {timeout:.0f}s",
                           "seconds": timeout})
                    s.kill(); s.start(); s.item = None
                if s.item is None and queue:
                    s.assign(queue.pop())
    for s in slots:
        try:
            s.parent.send(None)
        except Exception:  # noqa: BLE001
            pass
        s.kill()
    return {"done": n_done, "skipped": len(done_ids), "seconds": round(time.time() - t_start, 1)}
