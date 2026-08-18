#!/usr/bin/env python3
"""Phase 5 -- per-tissue transcriptional model.

Builds, per tissue: the expressed gene set, a tau tissue-specificity index, tissue-enriched
and tissue-specific classifications, and marker tables joined to the Swiss-Prot functional
annotation.

Operates on the rRNA-removed VST counts from 12_dge.R. Samples flagged as low-yield are
reported but, per the user's instruction, retained -- with their contribution to each tissue
mean made explicit so a tissue carried by one good replicate is visible as such.
"""
import argparse
import collections
import csv
import pathlib
import re
import statistics
import sys

ROOT = pathlib.Path(__file__).resolve().parents[1]


def read_matrix(path):
    with open(path) as fh:
        r = csv.reader(fh)
        header = next(r)
        samples = header[1:]
        genes, mat = [], []
        for row in r:
            genes.append(row[0])
            mat.append([float(x) for x in row[1:]])
    return samples, genes, mat


def tau(values):
    """Yanai tau: 0 = uniform across tissues, 1 = exclusive to one tissue."""
    mx = max(values)
    if mx <= 0:
        return float("nan")
    return sum(1 - (v / mx) for v in values) / (len(values) - 1)


def main():
    ap = argparse.ArgumentParser()
    # Linear normalised counts, NOT VST. tau is defined on linear expression; on VST values
    # the log compression pushes every value/max ratio towards 1 and tau collapses to ~0
    # (it returned 0 tissue-specific genes when first run on VST).
    ap.add_argument("--vst", default=str(ROOT / "results/dge_filtered/normalized_counts_rRNArm.csv"))
    ap.add_argument("--annot", default=str(ROOT / "results/annotation/pc915_functional.tsv"))
    ap.add_argument("--outdir", default=str(ROOT / "results/tissue_models"))
    ap.add_argument("--gtf", default=str(ROOT / "refs/PC9.15/PC9.15_genomic.gtf"),
                    help="stock GTF carrying gene_id/protein_id in CDS records")
    ap.add_argument("--tau-specific", type=float, default=0.85)
    ap.add_argument("--top", type=int, default=50)
    args = ap.parse_args()

    vst = pathlib.Path(args.vst)
    if not vst.exists():
        sys.exit(f"missing {vst} -- run 12_dge.R first")

    samples, genes, mat = read_matrix(str(vst))

    tissue_of, yield_of = {}, {}
    with (ROOT / "metadata" / "runsheet.csv").open() as fh:
        for r in csv.DictReader(fh):
            tissue_of[r["sample_name"]] = r["Factor Value[Tissue]"]
    ybudget = ROOT / "results" / "read_budget.csv"
    if ybudget.exists():
        with ybudget.open() as fh:
            for r in csv.DictReader(fh):
                yield_of[r["sample"]] = int(r["mRNA"])

    groups = collections.defaultdict(list)
    for i, s in enumerate(samples):
        groups[tissue_of.get(s, "?")].append(i)
    tissues = sorted(groups)
    print(f"tissues: {tissues}")
    for t in tissues:
        wells = [samples[i].split("_sample_")[-1] for i in groups[t]]
        ys = [yield_of.get(samples[i], 0) for i in groups[t]]
        print(f"  {t:<20} n={len(groups[t])}  wells={','.join(wells)}  "
              f"mRNA counts={[f'{y:,}' for y in ys]}")

    # per-tissue mean expression
    means = []
    for row in mat:
        means.append([statistics.mean(row[i] for i in groups[t]) for t in tissues])

    outdir = pathlib.Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    # The count matrix is keyed on GTF gene_id (locus tags, PTI98_*) while the functional
    # annotation is keyed on protein accessions (KAJ*). The GTF CDS records carry both, so
    # build gene_id -> protein_id from there; without it the marker join silently returns
    # nothing.
    gene2prot = {}
    for line in open(args.gtf):
        if line.startswith("#"):
            continue
        p = line.split("\t")
        if len(p) < 9 or p[2] != "CDS":
            continue
        g = re.search(r'gene_id "([^"]+)"', p[8])
        pr = re.search(r'protein_id "([^"]+)"', p[8])
        if g and pr:
            gene2prot.setdefault(g.group(1), pr.group(1))

    by_prot = {}
    ap_path = pathlib.Path(args.annot)
    if ap_path.exists():
        with ap_path.open() as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                by_prot[r["protein_id"]] = r
    annot = {g: by_prot[p] for g, p in gene2prot.items() if p in by_prot}
    print(f"\ngene->protein map: {len(gene2prot):,}; genes with annotation: {len(annot):,}")

    # tau + classification
    rows = []
    for g, m in zip(genes, means):
        t = tau(m)
        top_i = max(range(len(m)), key=lambda i: m[i])
        rows.append({
            "gene": g, "tau": t, "top_tissue": tissues[top_i],
            **{f"mean_{tissues[i]}": round(m[i], 3) for i in range(len(tissues))},
        })

    with (outdir / "tau_specificity.csv").open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)

    valid = [r for r in rows if r["tau"] == r["tau"]]
    specific = [r for r in valid if r["tau"] >= args.tau_specific]
    print(f"\ngenes scored: {len(valid):,}")
    print(f"tissue-specific (tau >= {args.tau_specific}): {len(specific):,}")
    per_t = collections.Counter(r["top_tissue"] for r in specific)
    for t in tissues:
        print(f"  {t:<20} {per_t.get(t, 0):>5}")

    # markers per tissue
    for t in tissues:
        cand = sorted((r for r in valid if r["top_tissue"] == t),
                      key=lambda r: -r["tau"])[:args.top]
        path = outdir / f"markers_{t.replace(' ', '_')}.csv"
        with path.open("w", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(["gene", "tau", f"mean_{t}", "protein_name", "EC", "KEGG"])
            for r in cand:
                a = annot.get(r["gene"], {})
                w.writerow([r["gene"], round(r["tau"], 3), r[f"mean_{t}"],
                            a.get("protein_name", ""), a.get("EC", ""), a.get("KEGG", "")])
        print(f"wrote {path}")

    print(f"\nwrote {outdir / 'tau_specificity.csv'}")


if __name__ == "__main__":
    sys.exit(main())
