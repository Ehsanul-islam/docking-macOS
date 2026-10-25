"""Typed errors (plan §4.4). Codes live in the JSON-RPC server-error range."""


class HelperError(Exception):
    code = -32000
    kind = "ToolFailure"

    def __init__(self, message: str, log: str = ""):
        super().__init__(message)
        self.message = message
        self.log = log

    def to_rpc(self) -> dict:
        return {
            "code": self.code,
            "message": self.message,
            "data": {"type": self.kind, "log": self.log},
        }


class InvalidInput(HelperError):
    code = -32001
    kind = "InvalidInput"


class UnsupportedFeature(HelperError):
    code = -32002
    kind = "UnsupportedFeature"


class ToolFailure(HelperError):
    code = -32003
    kind = "ToolFailure"


class Cancelled(HelperError):
    code = -32004
    kind = "Cancelled"
