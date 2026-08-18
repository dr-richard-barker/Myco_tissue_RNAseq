#!/usr/bin/env python3
"""Phase 1 -- summarise the reference test-map and apply the decision gate.

Parses HISAT2 --new-summary files written by 02_testmap.sh. Reports per reference the
mean unique ("aligned exactly 1 time") and multi ("aligned >1 times") rates. The plan's
gate is that the winning reference must exceed 60% uniquely mapped; below that no public
assembly represents this strain and de novo assembly becomes the fallback.
"""
import csv
import pathlib
import re
import statistics
import sys

ROOT = pathlib.Path(__file__).resolve().parents[1]
TESTMAP = ROOT / "qc" / "testmap"
GATE = 60.0

PATTERNS = {
    "total": r"Total reads:\s*(\d+)",
    "unal": r"Aligned 0 time:\s*(\d+)",
    "uniq": r"Aligned 1 time:\s*(\d+)",
    "multi": r"Aligned >1 times:\s*(\d+)",
}

tissue = {}
rs = ROOT / "metadata" / "runsheet.csv"
if rs.exists():
    with rs.open() as fh:
        for r in csv.DictReader(fh):
            tissue[r["sample_name"]] = r["Factor Value[Tissue]"]


def parse(path):
    text = path.read_text()
    vals = {}
    for key, pat in PATTERNS.items():
        m = re.search(pat, text)
        vals[key] = int(m.group(1)) if m else 0
    t = vals["total"] or 1
    return {
        "unique_pct": 100 * vals["uniq"] / t,
        "multi_pct": 100 * vals["multi"] / t,
        "unal_pct": 100 * vals["unal"] / t,
        "total": vals["total"],
    }


files = sorted(TESTMAP.glob("*.summary.txt"))
if not files:
    sys.exit(f"no HISAT2 summaries in {TESTMAP} -- has 02_testmap.sh finished?")

data = {}
for p in files:
    label, sample = p.name.split(".summary.txt")[0].split("__")
    data.setdefault(label, {})[sample] = parse(p)

print(f"{'reference':<10}{'n':>4}{'unique%':>10}{'multi%':>9}{'unaligned%':>12}{'total_aln%':>12}")
print("-" * 57)
summary = {}
for label, samples in sorted(
    data.items(), key=lambda kv: -statistics.mean(v["unique_pct"] for v in kv[1].values())
):
    uniq = [v["unique_pct"] for v in samples.values()]
    multi = [v["multi_pct"] for v in samples.values()]
    unal = [v["unal_pct"] for v in samples.values()]
    summary[label] = statistics.mean(uniq)
    print(f"{label:<10}{len(uniq):>4}{statistics.mean(uniq):>10.2f}{statistics.mean(multi):>9.2f}"
          f"{statistics.mean(unal):>12.2f}{100 - statistics.mean(unal):>12.2f}")

best = max(summary, key=summary.get)
runner_up = sorted(summary.values(), reverse=True)[1] if len(summary) > 1 else 0
print(f"\nbest reference: {best} ({summary[best]:.2f}% unique, "
      f"+{summary[best] - runner_up:.2f} pts over runner-up)")

if tissue:
    print(f"\nper-sample unique% for {best}:")
    by_tissue = {}
    for sample, v in data[best].items():
        by_tissue.setdefault(tissue.get(sample, "?"), []).append(
            (sample.split("_sample_")[-1], v["unique_pct"]))
    for t in sorted(by_tissue):
        vals = sorted(by_tissue[t])
        detail = "  ".join(f"{w}:{u:.1f}" for w, u in vals)
        print(f"  {t:<20} mean={statistics.mean([u for _, u in vals]):5.2f}   {detail}")

print()
if summary[best] < GATE:
    print(f"GATE FAILED: {summary[best]:.2f}% < {GATE}% unique. No public assembly represents "
          f"this strain well; de novo transcriptome assembly is the fallback (plan: Risks).")
    sys.exit(1)
print(f"GATE PASSED: {summary[best]:.2f}% >= {GATE}% unique. Proceed with {best}.")
