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


class NeedsTemplate(InvalidInput):
    """Bond orders cannot be established. The caller must supply a SMILES (or net charge if the file has
    hydrogens). `needs` tells the UI what to ask for."""
    code = -32005
    kind = "NeedsTemplate"

    def __init__(self, message: str, resname: str = "", needs: str = "smiles", log: str = ""):
        super().__init__(message, log)
        self.resname, self.needs = resname, needs

    def to_rpc(self) -> dict:
        d = super().to_rpc()
        d["data"].update({"resname": self.resname, "needs": self.needs})
        return d


class SelectionError(InvalidInput):
    """A selection expression could not be parsed or evaluated. `pos` is the character offset."""
    code = -32006
    kind = "SelectionError"

    def __init__(self, message: str, pos: int = 0, expr: str = ""):
        super().__init__(message)
        self.pos, self.expr = pos, expr

    def to_rpc(self) -> dict:
        d = super().to_rpc()
        d["data"].update({"pos": self.pos, "expr": self.expr})
        return d
