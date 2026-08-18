#!/usr/bin/env python3
"""Phase 6d -- reaction-level expression contrast across the metabolic network.

Motivation, and a correction to the obvious approach. Reaction presence/absence in RIPTiDe
context-specific models looked like a rich source of tissue-specific biology: 39 reactions
appear in the exudophore model and in none of the other three. But checking those reactions
against the expression that supposedly justifies them shows that most are not
exudophore-elevated at all -- several are expressed as highly or more highly in nodule, and 28
of the 39 have no gene association whatsoever, so their presence reflects network feasibility
rather than transcription.

Reaction inclusion in a context-specific model is therefore not, by itself, evidence of
tissue-specific metabolism. This script computes the evidence directly: map each model
reaction to its genes through the GPR, summarise expression per tissue, and rank reactions by
the contrast between the exudophore and every other tissue. Only reactions whose expression
actually separates survive.
"""
import argparse, collections, csv, pathlib, re, statistics, sys, warnings
warnings.filterwarnings("ignore")
import cobra

ROOT = pathlib.Path(__file__).resolve().parents[1]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--model", default=str(ROOT / "models/BOM_ss5_gapfilled.xml"))
    ap.add_argument("--expr", default=str(ROOT / "results/normalised_all_genes.csv"))
    ap.add_argument("--gtf", default=str(ROOT / "refs/BOM_ss5/BOM_ss5_genomic.gtf"))
    ap.add_argument("--out", default=str(ROOT / "results/tissue_metabolism/reaction_expression.csv"))
    args = ap.parse_args()

    p2g = {}
    for line in open(args.gtf):
        if line.startswith("#"): continue
        p = line.split("\t")
        if len(p) < 9 or p[2] != "CDS": continue
        g = re.search(r'gene_id "([^"]+)"', p[8]); pr = re.search(r'protein_id "([^"]+)"', p[8])
        if g and pr: p2g.setdefault(pr.group(1), g.group(1))

    with open(args.expr) as fh:
        r = csv.reader(fh); hdr = next(r); samples = hdr[1:]
        expr = {row[0]: [float(x) for x in row[1:]] for row in r}
    tis = {x["sample_name"]: x["Factor Value[Tissue]"]
           for x in csv.DictReader((ROOT / "metadata/runsheet.csv").open())}
    grp = collections.defaultdict(list)
    for i, s in enumerate(samples): grp[tis[s]].append(i)
    tissues = sorted(grp)

    model = cobra.io.read_sbml_model(args.model)
    rows = []
    for rx in model.reactions:
        accs = re.findall(r"KAN[0-9.]+", rx.gene_reaction_rule or "")
        genes = [p2g.get(a) for a in accs]
        genes = [g for g in genes if g in expr]
        if not genes: continue
        # reaction expression = max over isoenzymes (any one suffices to catalyse)
        per = {}
        for t in tissues:
            per[t] = max(statistics.mean(expr[g][i] for i in grp[t]) for g in genes)
        others = [per[t] for t in tissues if t != "Exudophore"]
        exu = per["Exudophore"]
        rows.append(dict(reaction=rx.id, name=(rx.name or rx.id)[:70], n_genes=len(genes),
                         **{t.replace(" ", "_"): round(per[t], 2) for t in tissues},
                         max_other=round(max(others), 2),
                         ratio=round(exu / max(max(others), 1e-9), 2)))

    rows.sort(key=lambda r: -r["ratio"])
    with open(args.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)

    print(f"model reactions with gene evidence and expression: {len(rows)}\n")
    print(f"{'reaction':<11}{'enzyme':<46}{'Exudo':>8}{'maxOther':>10}{'ratio':>7}")
    print("-" * 84)
    for r in rows[:18]:
        print(f"{r['reaction']:<11}{r['name'][:44]:<46}{r['Exudophore']:>8.1f}"
              f"{r['max_other']:>10.1f}{r['ratio']:>7.1f}")
    strong = [r for r in rows if r["ratio"] >= 2 and r["Exudophore"] >= 5]
    print(f"\nreactions >=2x the highest other tissue and above a floor of 5: {len(strong)}")
    print(f"wrote {args.out}")


if __name__ == "__main__":
    sys.exit(main())
