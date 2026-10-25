"""JSON-RPC 2.0 server over stdin/stdout (plan §4.4).

One request per line. Responses carry result, warnings, provenance and
artifacts. stdout is reserved for protocol traffic; logs go to stderr.
"""
import json
import sys
import threading

from . import __version__
from .errors import HelperError, InvalidInput, ToolFailure
from .methods import METHODS


class Server:
    def __init__(self, stdin=None, stdout=None):
        self.stdin = stdin or sys.stdin
        self.stdout = stdout or sys.stdout
        self._lock = threading.Lock()
        self.cancelled: set = set()

    def send(self, msg: dict) -> None:
        with self._lock:
            self.stdout.write(json.dumps(msg) + "\n")
            self.stdout.flush()

    def progress(self, job_id, fraction: float, message: str = "") -> None:
        self.send({"jsonrpc": "2.0", "method": "progress",
                   "params": {"job_id": job_id, "fraction": fraction, "message": message}})

    def handle(self, req: dict) -> dict | None:
        rid = req.get("id")
        method = req.get("method")
        params = req.get("params") or {}
        if method == "cancel":  # notification: mark job cancelled
            self.cancelled.add(params.get("job_id"))
            return None
        try:
            if req.get("jsonrpc") != "2.0" or not isinstance(method, str):
                raise InvalidInput("malformed JSON-RPC request")
            fn = METHODS.get(method)
            if fn is None:
                return {"jsonrpc": "2.0", "id": rid,
                        "error": {"code": -32601, "message": f"unknown method {method}"}}
            out = fn(params, self)
            result = {"result": None, "warnings": [], "provenance": {}, "artifacts": []}
            result.update(out)
            return {"jsonrpc": "2.0", "id": rid, "result": result}
        except HelperError as e:
            return {"jsonrpc": "2.0", "id": rid, "error": e.to_rpc()}
        except Exception as e:  # never let a tool crash the loop
            return {"jsonrpc": "2.0", "id": rid, "error": ToolFailure(repr(e)).to_rpc()}

    def serve(self) -> None:
        for line in self.stdin:
            line = line.strip()
            if not line:
                continue
            try:
                req = json.loads(line)
            except json.JSONDecodeError as e:
                self.send({"jsonrpc": "2.0", "id": None,
                           "error": {"code": -32700, "message": f"parse error: {e}"}})
                continue
            resp = self.handle(req)
            if resp is not None and "id" in req:
                self.send(resp)


def main() -> None:
    print(f"dockchem-helper {__version__} ready", file=sys.stderr)
    Server().serve()


if __name__ == "__main__":
    main()
