# Phase 0 status

Exit criteria (plan §18): §12.1 suite passes on both architectures; helper runs RPC round-trip.

| Item | State |
|---|---|
| Helper RPC round-trip, typed errors, progress, cancel | Done, tested (pytest) |
| `.dockproj` store + provenance log + checksum verify | Done, tested |
| Equivalence harness | Done, tested against stub engines only |
| Vina build script | Written, **unrun** (needs macOS + per-arch static Boost) |
| Swift package | Written, **uncompiled** (needs macOS: `swift test --package-path swift`) |
| Python packaging/signing prototype (D-2) | Not started |
| 50-complex equivalence data set | Not started |
| Cross-architecture numeric measurement (R1) | Not started |

Open [VERIFY]s touched: upstream Vina CLI builds via Makefile (not CMake); Boost static per-arch.
