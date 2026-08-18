#!/usr/bin/env python3
"""Phase 3/4 -- read-budget and per-sample QC from the full-depth counts.

Produces the table that decides which samples can carry a differential-expression result.
Splits assigned counts by gene_biotype so the rRNA burden is explicit, and reports genes
detected at two thresholds as a direct proxy for usable complexity.
"""
import argparse
import collections
import csv
import json
import pathlib
import re
import statistics
import sys

ROOT = pathlib.Path(__file__).resolve().parents[1]


def biotypes(gtf):
    bt = {}
    for line in open(gtf):
        if line.startswith("#"):
            continue
        p = line.split("\t")
        if len(p) < 9 or p[2] != "gene":
            continue
        g = re.search(r'gene_id "([^"]+)"', p[8])
        b = re.search(r'gene_biotype "([^"]+)"', p[8])
        if g:
            bt[g.group(1)] = b.group(1) if b else "unknown"
    return bt


def load_counts(path):
    samples, rows = None, []
    for line in open(path):
        if line.startswith("#"):
            continue
        p = line.rstrip("\n").split("\t")
        if p[0] == "Geneid":
            samples = [pathlib.Path(x).name.split(".")[0] for x in p[6:]]
            continue
        rows.append((p[0], [int(x) for x in p[6:]]))
    return samples, rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--counts", default=str(ROOT / "counts" / "counts_dedup.txt"))
    ap.add_argument("--gtf", default=str(ROOT / "refs" / "PC9.15" / "PC9.15_final.gtf"))
    ap.add_argument("--min-mrna", type=int, default=100_000,
                    help="assigned mRNA counts below which a sample is flagged")
    args = ap.parse_args()

    bt = biotypes(args.gtf)
    samples, rows = load_counts(args.counts)
    if not samples:
        sys.exit(f"no samples parsed from {args.counts}")

    tissue, raw_reads = {}, {}
    with (ROOT / "metadata" / "runsheet.csv").open() as fh:
        for r in csv.DictReader(fh):
            tissue[r["sample_name"]] = r["Factor Value[Tissue]"]

    for s in samples:
        j = ROOT / "qc" / "fastp" / f"{s}.json"
        if j.exists():
            raw_reads[s] = json.load(j.open())["summary"]["before_filtering"]["total_reads"]

    agg = collections.defaultdict(lambda: [0] * len(samples))
    for g, v in rows:
        k = bt.get(g, "unknown")
        k = "rRNA" if k == "rRNA" else ("tRNA" if k == "tRNA" else "mRNA")
        for i, x in enumerate(v):
            agg[k][i] += x

    out = []
    for i, s in enumerate(samples):
        mrna, rrna = agg["mRNA"][i], agg["rRNA"][i]
        det1 = sum(1 for g, v in rows if bt.get(g) == "protein_coding" and v[i] > 0)
        det10 = sum(1 for g, v in rows if bt.get(g) == "protein_coding" and v[i] >= 10)
        out.append({
            "sample": s, "well": s.split("_sample_")[-1],
            "tissue": tissue.get(s, "?"), "raw": raw_reads.get(s, 0),
            "rRNA": rrna, "mRNA": mrna, "det1": det1, "det10": det10,
            "pct": 100 * mrna / raw_reads[s] if raw_reads.get(s) else float("nan"),
        })

    out.sort(key=lambda r: (r["tissue"], r["well"]))
    print(f"{'tissue':<19}{'well':<5}{'raw reads':>12}{'rRNA':>12}{'mRNA':>12}"
          f"{'mRNA%raw':>10}{'genes>=1':>10}{'>=10':>8}  flag")
    print("-" * 98)
    for r in out:
        flag = "LOW" if r["mRNA"] < args.min_mrna else ""
        print(f"{r['tissue']:<19}{r['well']:<5}{r['raw']:>12,}{r['rRNA']:>12,}{r['mRNA']:>12,}"
              f"{r['pct']:>9.1f}%{r['det1']:>10,}{r['det10']:>8,}  {flag}")

    print("-" * 98)
    by_t = collections.defaultdict(list)
    for r in out:
        by_t[r["tissue"]].append(r)
    for t, rs in sorted(by_t.items()):
        usable = sum(1 for r in rs if r["mRNA"] >= args.min_mrna)
        print(f"{t:<19}{'mean':<5}{statistics.mean(r['raw'] for r in rs):>12,.0f}"
              f"{statistics.mean(r['rRNA'] for r in rs):>12,.0f}"
              f"{statistics.mean(r['mRNA'] for r in rs):>12,.0f}"
              f"{statistics.mean(r['pct'] for r in rs):>9.1f}%"
              f"{statistics.mean(r['det1'] for r in rs):>10,.0f}"
              f"{statistics.mean(r['det10'] for r in rs):>8,.0f}"
              f"  n usable = {usable}/{len(rs)}")

    dest = ROOT / "results" / "read_budget.csv"
    dest.parent.mkdir(exist_ok=True)
    with dest.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(out[0]))
        w.writeheader()
        w.writerows(out)
    print(f"\nwrote {dest}")


if __name__ == "__main__":
    sys.exit(main())
