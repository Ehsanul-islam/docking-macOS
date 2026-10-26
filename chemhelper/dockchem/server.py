"""JSON-RPC 2.0 server over stdin/stdout (plan §4.4).

One request per line. Responses carry result, warnings, provenance and artifacts. stdout is reserved for protocol
traffic; logs go to stderr.

Requests run one at a time on a worker thread, in order; the reading thread stays free, so a `cancel` notification is
seen while a job is running and trips that job's cancellation token (see cancel.py). A cancel that arrives before its
job starts cancels it on arrival.
"""
import json
import queue
import sys
import threading
import traceback

from . import __version__, cancel
from .errors import HelperError, InternalError, InvalidInput
from .methods import METHODS


class Server:
    def __init__(self, stdin=None, stdout=None):
        self.stdin = stdin or sys.stdin
        self.stdout = stdout or sys.stdout
        self._lock = threading.Lock()
        self.cancelled: set = set()
        self.tokens: dict = {}

    def send(self, msg: dict) -> None:
        with self._lock:
            self.stdout.write(json.dumps(msg) + "\n")
            self.stdout.flush()

    def progress(self, job_id, fraction: float, message: str = "") -> None:
        self.send({"jsonrpc": "2.0", "method": "progress",
                   "params": {"job_id": job_id, "fraction": fraction, "message": message}})

    def _cancel(self, job_id) -> None:
        self.cancelled.add(job_id)
        tok = self.tokens.get(job_id)
        if tok is not None:
            tok.cancel()

    def handle(self, req: dict) -> dict | None:
        rid = req.get("id") if isinstance(req, dict) else None
        method = req.get("method") if isinstance(req, dict) else None
        params = (req.get("params") if isinstance(req, dict) else None) or {}
        if not isinstance(params, dict):
            return {"jsonrpc": "2.0", "id": rid, "error": InvalidInput("params must be an object").to_rpc()}
        if method == "cancel":
            self._cancel(params.get("job_id"))
            return None
        job_id = params.get("job_id", rid)
        token = cancel.Token()
        if job_id in self.cancelled:
            token.cancel()
        self.tokens[job_id] = token
        handle = cancel.bind(token)
        try:
            if not isinstance(req, dict) or req.get("jsonrpc") != "2.0" or not isinstance(method, str):
                raise InvalidInput("malformed JSON-RPC request")
            fn = METHODS.get(method)
            if fn is None:
                return {"jsonrpc": "2.0", "id": rid,
                        "error": {"code": -32601, "message": f"unknown method {method}"}}
            token.raise_if_cancelled()
            out = fn(params, self)
            result = {"result": None, "warnings": [], "provenance": {}, "artifacts": []}
            result.update(out)
            return {"jsonrpc": "2.0", "id": rid, "result": result}
        except HelperError as e:
            return {"jsonrpc": "2.0", "id": rid, "error": e.to_rpc()}
        except Exception as e:  # never let a tool crash the loop; but say plainly that this is a bug
            tb = traceback.format_exc()
            print(tb, file=sys.stderr)
            return {"jsonrpc": "2.0", "id": rid, "error": InternalError(
                f"internal error in {method}: {type(e).__name__}: {e}", log=tb[-2000:]).to_rpc()}
        finally:
            cancel.unbind(handle)
            self.tokens.pop(job_id, None)

    def _worker(self, jobs: "queue.Queue") -> None:
        while True:
            req = jobs.get()
            if req is None:
                return
            if isinstance(req, dict) and "__immediate__" in req:       # e.g. a parse error: sent in request order
                self.send(req["__immediate__"])
                continue
            resp = self.handle(req)
            if resp is not None and isinstance(req, dict) and "id" in req:
                self.send(resp)

    def serve(self) -> None:
        jobs: queue.Queue = queue.Queue()
        worker = threading.Thread(target=self._worker, args=(jobs,), daemon=True)
        worker.start()
        for line in self.stdin:
            line = line.strip()
            if not line:
                continue
            try:
                req = json.loads(line)
            except json.JSONDecodeError as e:
                jobs.put({"__immediate__": {"jsonrpc": "2.0", "id": None,
                                            "error": {"code": -32700, "message": f"parse error: {e}"}}})
                continue
            if isinstance(req, dict) and req.get("method") == "cancel":     # act at once, even mid-job
                self._cancel((req.get("params") or {}).get("job_id") if isinstance(req.get("params"), dict) else None)
                continue
            jobs.put(req)
        jobs.put(None)
        worker.join()


def main() -> None:
    print(f"dockchem-helper {__version__} ready", file=sys.stderr)
    Server().serve()


if __name__ == "__main__":
    main()
