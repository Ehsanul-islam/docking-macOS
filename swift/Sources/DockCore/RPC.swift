import Foundation

/// JSON-RPC 2.0 wire types for the chemistry helper (plan §4.4).
public struct RPCRequest: Encodable, Sendable {
    public let jsonrpc = "2.0"
    public let id: Int
    public let method: String
    public let params: [String: JSONValue]
    public init(id: Int, method: String, params: [String: JSONValue] = [:]) {
        self.id = id; self.method = method; self.params = params
    }
}

/// The helper's typed error kinds. Unknown kinds (newer helper) are kept as raw text, never a decoding failure.
public struct HelperErrorKind: RawRepresentable, Hashable, Decodable, Sendable {
    public let rawValue: String
    public init(rawValue: String) { self.rawValue = rawValue }
    public init(from decoder: Decoder) throws { rawValue = try decoder.singleValueContainer().decode(String.self) }

    public static let invalidInput = HelperErrorKind(rawValue: "InvalidInput")
    public static let unsupportedFeature = HelperErrorKind(rawValue: "UnsupportedFeature")
    public static let toolFailure = HelperErrorKind(rawValue: "ToolFailure")
    public static let cancelled = HelperErrorKind(rawValue: "Cancelled")
    public static let needsTemplate = HelperErrorKind(rawValue: "NeedsTemplate")
    public static let selectionError = HelperErrorKind(rawValue: "SelectionError")
    public static let unusableStructure = HelperErrorKind(rawValue: "UnusableStructure")
    public static let formatLimit = HelperErrorKind(rawValue: "FormatLimit")
    public static let internalError = HelperErrorKind(rawValue: "InternalError")
}

public struct RPCError: Decodable, Error, Sendable {
    public struct Details: Decodable, Sendable {
        public let type: HelperErrorKind?
        public let log: String?
        /// Extra, kind-specific fields (e.g. `assessment` for UnusableStructure, `needs` for NeedsTemplate).
        public let extra: [String: JSONValue]

        public init(from decoder: Decoder) throws {
            let all = try decoder.singleValueContainer().decode([String: JSONValue].self)
            type = all["type"]?.string.map { HelperErrorKind(rawValue: $0) }
            log = all["log"]?.string
            extra = all.filter { $0.key != "type" && $0.key != "log" }
        }
    }
    public let code: Int
    public let message: String
    public let data: Details?

    public var kind: HelperErrorKind? { data?.type }
    public var isCancelled: Bool { kind == .cancelled }
}

extension RPCError: LocalizedError {
    public var errorDescription: String? { message }
}

public struct RPCResponse: Decodable, Sendable {
    public let id: Int?
    public let error: RPCError?
    public let result: JSONValue?
}

/// What a successful call returns: the helper's envelope (plan §4.4).
public struct HelperResult: Sendable {
    public let result: JSONValue
    public let warnings: [String]
    public let provenance: JSONValue
    public let artifacts: [String]
}
