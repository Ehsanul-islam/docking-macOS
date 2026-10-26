#!/usr/bin/env python3
"""Summarise a protein_batch run into notes: what happened, where, and (from the recorded evidence) why."""
import collections
import json
import re
import sys

rows = [json.loads(l) for l in open(sys.argv[1] + "/proteins.jsonl")]
n = len(rows)
norm = lambda s: re.sub(r"\d+", "N", (s or "").replace("\n", " "))[:110]  # noqa: E731
pct = lambda a, b: f"{a}/{b} = {100 * a / max(1, b):.1f}%"  # noqa: E731
print(f"# Survey notes ({n} entries)\n")
print("## Outcome overall")
for c, k in collections.Counter(r["category"] for r in rows).most_common():
    print(f"- {c}: {pct(k, n)}")
print("\n## Outcome by organism group and method (ok / degraded / unusable / failed+exception / timeout)")
for key in ("organism_group", "method"):
    grp = collections.defaultdict(collections.Counter)
    for r in rows:
        c = r["category"]
        grp[r.get(key)][("failed" if c in ("failed", "exception", "crash") else c)] += 1
    for g, cnt in sorted(grp.items(), key=lambda kv: str(kv[0])):
        t = sum(cnt.values())
        print(f"- {key}={g}: n={t} " + ", ".join(f"{c} {cnt[c]}" for c in ("ok", "degraded", "unusable", "failed", "timeout") if cnt[c]))
print("\n## Assessment verdict vs outcome")
va = collections.defaultdict(collections.Counter)
for r in rows:
    va[r.get("verdict", "none")][r["category"]] += 1
for v, cnt in va.items():
    print(f"- {v}: {dict(cnt)}")
print("\n## Failures, grouped by message (where/why)")
g = collections.defaultdict(list)
for r in rows:
    if r["category"] in ("failed", "exception", "crash", "missing_file"):
        g[(r["category"], r.get("error_kind", ""), norm(r.get("error")))].append(r)
for (c, kind, e), rs in sorted(g.items(), key=lambda kv: -len(kv[1])):
    print(f"- [{len(rs)}] {c} {kind}: {e}\n    ids: {', '.join(x['id'] for x in rs[:8])}  methods: {dict(collections.Counter(x.get('method') for x in rs))}")
    tail = (rs[0].get("trace") or rs[0].get("log_tail") or "").strip().splitlines()
    if tail:
        print("    evidence: " + " | ".join(l.strip()[:120] for l in tail[-3:]))
print("\n## Timeouts (15 min limit)")
for r in rows:
    if r["category"] == "timeout":
        print(f"- {r['id']} method={r.get('method')} group={r.get('organism_group')} atoms={r.get('n_atoms')} verdict={r.get('verdict')}")
print("\n## Refused / unusable reasons")
u = collections.defaultdict(list)
for r in rows:
    if r["category"] == "unusable":
        u[norm(r.get("error"))[:100]].append(r["id"])
for e, ids in sorted(u.items(), key=lambda kv: -len(kv[1])):
    print(f"- [{len(ids)}] {e}  e.g. {', '.join(ids[:6])}")
print("\n## Degraded: what was left out, how far from the site")
d = collections.Counter()
for r in rows:
    if r["category"] == "degraded":
        for x in r.get("dropped", []):
            d[norm(x.split(" :: ")[1])[:80]] += 1
for e, k in d.most_common(12):
    print(f"- [{k}] {e}")
print("\n## Prepared entries: writer and cross-check")
prep = [r for r in rows if r["category"] in ("ok", "degraded")]
print("- prepared:", pct(len(prep), n), "| writer:", dict(collections.Counter(r.get("writer") for r in prep)),
      "| cross-check:", dict(collections.Counter(r.get("crosscheck") for r in prep)))
for r in prep:
    if r.get("crosscheck") == "disagree":
        print(f"  - disagree {r['id']}: {r['crosscheck_examples'][:1]}")
    if r.get("crosscheck") == "unavailable":
        pass
un = collections.Counter(norm(r.get("crosscheck_reason")) for r in prep if r.get("crosscheck") == "unavailable")
if un:
    print("- Meeko cross-check unavailable:", dict(un.most_common(4)))
print("\n## Most frequent warnings (prepared entries)")
w = collections.Counter()
for r in prep:
    for x in set(norm(y)[:90] for y in r.get("warnings", [])):
        w[x] += 1
for e, k in w.most_common(14):
    print(f"- [{k}] {e}")
print("\n## Assessment findings (fix/block) frequency")
f = collections.Counter(c for r in rows for c in r.get("assessment_findings", []))
print("- " + ", ".join(f"{c} {k}" for c, k in f.most_common()))
print("\n## Time")
t = sorted(r.get("seconds", 0) for r in rows if r["category"] != "timeout")
print(f"- median {t[len(t) // 2]:.0f}s, p90 {t[int(len(t) * .9)]:.0f}s, p99 {t[int(len(t) * .99)]:.0f}s, max {t[-1]:.0f}s")
slow = sorted(rows, key=lambda r: -r.get("seconds", 0))[:8]
print("- slowest: " + ", ".join(f"{r['id']}({r.get('n_atoms')} atoms, {r.get('seconds', 0):.0f}s, {r['category']})" for r in slow))
