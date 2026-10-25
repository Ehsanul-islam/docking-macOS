import Foundation

/// JSON-RPC 2.0 wire types for the chemistry helper (plan §4.4).
public struct RPCRequest: Encodable {
    public let jsonrpc = "2.0"
    public let id: Int
    public let method: String
    public let params: [String: String]
    public init(id: Int, method: String, params: [String: String] = [:]) {
        self.id = id; self.method = method; self.params = params
    }
}

public enum HelperErrorKind: String, Decodable {
    case invalidInput = "InvalidInput"
    case unsupportedFeature = "UnsupportedFeature"
    case toolFailure = "ToolFailure"
    case cancelled = "Cancelled"
}

public struct RPCError: Decodable, Error {
    public struct Data: Decodable { public let type: HelperErrorKind?; public let log: String? }
    public let code: Int
    public let message: String
    public let data: Data?
}

public struct RPCResponse: Decodable {
    public let id: Int?
    public let error: RPCError?
}
