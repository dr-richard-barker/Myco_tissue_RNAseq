#!/usr/bin/env python3
"""Phase 5 -- test whether tau markers are group-wide or single-replicate artefacts.

tau is computed on tissue MEANS, so one outlying replicate can manufacture a high-tau
"marker" for a whole tissue. This applies the stricter test the mean cannot: a marker counts
as supported only if EVERY replicate of its tissue exceeds the highest value seen in any
sample of any other tissue. Reported per tissue as supported/total.
"""
import csv, sys, pathlib

ROOT = pathlib.Path(__file__).resolve().parents[1]
norm = sys.argv[1] if len(sys.argv) > 1 else str(ROOT / "results/dge_BOM_ss5_filtered/normalized_counts_rRNArm.csv")
mdir = pathlib.Path(sys.argv[2] if len(sys.argv) > 2 else ROOT / "results/tissue_models_BOM_ss5")
topn = int(sys.argv[3]) if len(sys.argv) > 3 else 50

M = {}
with open(norm) as fh:
    r = csv.reader(fh); hdr = next(r); samples = hdr[1:]
    for row in r: M[row[0]] = [float(x) for x in row[1:]]

tis = {x["sample_name"]: x["Factor Value[Tissue]"] for x in csv.DictReader(open(ROOT / "metadata/runsheet.csv"))}
groups = {}
for i, s in enumerate(samples): groups.setdefault(tis[s], []).append(i)

out = []
for t, idx in sorted(groups.items()):
    f = mdir / f"markers_{t.replace(' ', '_')}.csv"
    if not f.exists(): continue
    rows = list(csv.DictReader(open(f)))[:topn]
    ok = tot = 0
    for x in rows:
        v = M.get(x["gene"])
        if not v: continue
        tot += 1
        others = [v[i] for j, g in groups.items() if j != t for i in g]
        if min(v[i] for i in idx) > (max(others) if others else 0):
            ok += 1
            out.append((t, x["gene"], x.get("protein_name", ""), x["tau"]))
    print(f"{t:<20} n={len(idx)}  {ok:>3}/{tot:<3} markers supported by ALL replicates "
          f"({100*ok/max(tot,1):.0f}%)")

dest = mdir / "markers_robust.csv"
with dest.open("w", newline="") as fh:
    w = csv.writer(fh); w.writerow(["tissue", "gene", "protein_name", "tau"]); w.writerows(out)
print(f"\nwrote {dest}  ({len(out)} robust markers)")
