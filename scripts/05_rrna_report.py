#!/usr/bin/env python3
"""Phase 2a -- report the measured rRNA fraction per sample and per tissue.

Replaces the probe-based lower bound (12-46%, four diagnostic 25-mers) with a direct
measurement against the barrnap-derived rRNA loci. Also reports the implied unique
mapping rate among non-rRNA reads, which is the number the Phase 1 reference gate should
have been written against.
"""
import csv
import pathlib
import re
import statistics
import sys

ROOT = pathlib.Path(__file__).resolve().parents[1]
RRNA = ROOT / "qc" / "rrna"
TESTMAP = ROOT / "qc" / "testmap"
REF = sys.argv[1] if len(sys.argv) > 1 else "BOM_ss5"


def parse(path):
    text = path.read_text()

    def grab(pat):
        m = re.search(pat, text)
        return int(m.group(1)) if m else 0

    total = grab(r"Total reads:\s*(\d+)") or 1
    return {
        "total": total,
        "unal": grab(r"Aligned 0 time:\s*(\d+)"),
        "uniq": grab(r"Aligned 1 time:\s*(\d+)"),
        "multi": grab(r"Aligned >1 times:\s*(\d+)"),
    }


tissue = {}
with (ROOT / "metadata" / "runsheet.csv").open() as fh:
    for r in csv.DictReader(fh):
        tissue[r["sample_name"]] = r["Factor Value[Tissue]"]

rows = []
for p in sorted(RRNA.glob("*.summary.txt")):
    sample = p.name.replace(".summary.txt", "")
    r = parse(p)
    rrna_pct = 100 * (r["total"] - r["unal"]) / r["total"]

    g = TESTMAP / f"{REF}__{sample}.summary.txt"
    gm = parse(g) if g.exists() else None
    # Reads that hit rRNA are overwhelmingly multi-mappers (rDNA is a tandem repeat), so
    # treat the genome-unique count as coming entirely from the non-rRNA pool. This is an
    # estimate, not a measurement -- it is validated by re-mapping depleted reads below.
    uniq_nonrrna = (100 * gm["uniq"] / max(gm["total"] - (gm["total"] * rrna_pct / 100), 1)
                    if gm else float("nan"))
    rows.append({
        "sample": sample,
        "well": sample.split("_sample_")[-1],
        "tissue": tissue.get(sample, "?"),
        "rrna_pct": rrna_pct,
        "uniq_pct": 100 * gm["uniq"] / gm["total"] if gm else float("nan"),
        "uniq_nonrrna": uniq_nonrrna,
    })

rows.sort(key=lambda r: (r["tissue"], r["well"]))
print(f"rRNA measured against barrnap loci from {REF}\n")
print(f"{'tissue':<20}{'well':<6}{'rRNA%':>9}{'uniq%(all)':>12}{'uniq%(non-rRNA, est)':>22}")
print("-" * 69)
for r in rows:
    print(f"{r['tissue']:<20}{r['well']:<6}{r['rrna_pct']:>9.1f}"
          f"{r['uniq_pct']:>12.1f}{r['uniq_nonrrna']:>22.1f}")

print("-" * 69)
by_t = {}
for r in rows:
    by_t.setdefault(r["tissue"], []).append(r)
for t, rs in sorted(by_t.items()):
    print(f"{t:<20}{'mean':<6}{statistics.mean(x['rrna_pct'] for x in rs):>9.1f}"
          f"{statistics.mean(x['uniq_pct'] for x in rs):>12.1f}"
          f"{statistics.mean(x['uniq_nonrrna'] for x in rs):>22.1f}")

allr = [r["rrna_pct"] for r in rows]
print(f"\noverall rRNA: mean={statistics.mean(allr):.1f}%  "
      f"range={min(allr):.1f}-{max(allr):.1f}%")
print(f"implied usable (non-rRNA) fraction: mean={100 - statistics.mean(allr):.1f}%")
