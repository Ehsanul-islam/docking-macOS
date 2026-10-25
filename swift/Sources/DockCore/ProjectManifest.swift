import Foundation

/// Mirrors schemas/project.schema.json (plan §4.6).
public struct ProjectManifest: Codable, Equatable {
    public struct Input: Codable, Equatable {
        public var path: String
        public var sha256: String
    }
    public var schemaVersion: Int
    public var id: String
    public var name: String
    public var inputs: [Input]

    enum CodingKeys: String, CodingKey {
        case schemaVersion = "schema_version", id, name, inputs
    }
}
