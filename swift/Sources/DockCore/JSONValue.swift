import Foundation

/// Any JSON value. Used for request parameters and for results whose shape the app reads field by field.
public enum JSONValue: Codable, Equatable, Sendable, ExpressibleByStringLiteral, ExpressibleByIntegerLiteral,
                       ExpressibleByFloatLiteral, ExpressibleByBooleanLiteral, ExpressibleByArrayLiteral,
                       ExpressibleByDictionaryLiteral, ExpressibleByNilLiteral {
    case null
    case bool(Bool)
    case int(Int)
    case double(Double)
    case string(String)
    case array([JSONValue])
    case object([String: JSONValue])

    public init(stringLiteral v: String) { self = .string(v) }
    public init(integerLiteral v: Int) { self = .int(v) }
    public init(floatLiteral v: Double) { self = .double(v) }
    public init(booleanLiteral v: Bool) { self = .bool(v) }
    public init(arrayLiteral elements: JSONValue...) { self = .array(elements) }
    public init(dictionaryLiteral elements: (String, JSONValue)...) {
        self = .object(Dictionary(elements, uniquingKeysWith: { $1 }))
    }
    public init(nilLiteral: ()) { self = .null }

    public init(from decoder: Decoder) throws {
        let c = try decoder.singleValueContainer()
        if c.decodeNil() { self = .null }
        else if let v = try? c.decode(Bool.self) { self = .bool(v) }
        else if let v = try? c.decode(Int.self) { self = .int(v) }
        else if let v = try? c.decode(Double.self) { self = .double(v) }
        else if let v = try? c.decode(String.self) { self = .string(v) }
        else if let v = try? c.decode([JSONValue].self) { self = .array(v) }
        else if let v = try? c.decode([String: JSONValue].self) { self = .object(v) }
        else { throw DecodingError.dataCorruptedError(in: c, debugDescription: "unsupported JSON value") }
    }

    public func encode(to encoder: Encoder) throws {
        var c = encoder.singleValueContainer()
        switch self {
        case .null: try c.encodeNil()
        case .bool(let v): try c.encode(v)
        case .int(let v): try c.encode(v)
        case .double(let v): try c.encode(v)
        case .string(let v): try c.encode(v)
        case .array(let v): try c.encode(v)
        case .object(let v): try c.encode(v)
        }
    }

    // MARK: reading

    public subscript(key: String) -> JSONValue? {
        if case .object(let o) = self { return o[key] }
        return nil
    }
    public subscript(index: Int) -> JSONValue? {
        if case .array(let a) = self, a.indices.contains(index) { return a[index] }
        return nil
    }
    public var string: String? { if case .string(let v) = self { return v }; return nil }
    public var bool: Bool? { if case .bool(let v) = self { return v }; return nil }
    public var int: Int? {
        switch self {
        case .int(let v): return v
        case .double(let v) where v == v.rounded() && abs(v) < 1e15: return Int(v)
        default: return nil
        }
    }
    public var double: Double? {
        switch self {
        case .double(let v): return v
        case .int(let v): return Double(v)
        default: return nil
        }
    }
    public var array: [JSONValue]? { if case .array(let v) = self { return v }; return nil }
    public var object: [String: JSONValue]? { if case .object(let v) = self { return v }; return nil }
    public var isNull: Bool { if case .null = self { return true }; return false }
}
