#!/usr/bin/env python3
"""Phase 5 -- GO (CC/MF/BP) and KEGG enrichment per WGCNA module, then a readable title.

Hypergeometric test of each module against the background of all network genes carrying at
least one annotation. FDR is applied with Benjamini-Hochberg WITHIN each module x ontology,
which is the relevant family: the question asked of each module is "which terms of this
ontology are over-represented in it".

Caveats carried into the output:
  * only ~41% of the proteome has a Swiss-Prot hit, so annotation coverage is the limiting
    factor and absence of enrichment is weak evidence;
  * WGCNA was run at n=12, below the ~15 its authors recommend, so modules are exploratory;
  * module titles are a readable summary of the top enriched terms, not a claim of function.
"""
import argparse, collections, csv, math, pathlib, sys
from scipy.stats import hypergeom
from statsmodels.stats.multitest import multipletests

ROOT = pathlib.Path(__file__).resolve().parents[1]
ONTOLOGIES = ["GO_BP", "GO_MF", "GO_CC", "KEGG"]

# Words that make a term uninformative as a title.
VAGUE = ("binding", "catalytic activity", "metabolic process", "biosynthetic process",
         "cellular process", "cytoplasm", "membrane", "nucleus", "protein binding",
         "ATP binding", "metal ion binding", "cytosol", "intracellular")


def load_terms(path):
    g2t, t2g, tname = collections.defaultdict(set), collections.defaultdict(set), {}
    with open(path) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            key = (r["ontology"], r["term_id"])
            g2t[r["gene_id"]].add(key)
            t2g[key].add(r["gene_id"])
            tname[key] = r["term_name"]
    return g2t, t2g, tname


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--modules", default=str(ROOT / "results/wgcna/gene_modules.csv"))
    ap.add_argument("--terms", default=str(ROOT / "results/annotation/bom_ss5_terms.tsv"))
    ap.add_argument("--cazy", default=str(ROOT / "results/annotation/bom_ss5_cazymes.tsv"))
    ap.add_argument("--secretome", default=str(ROOT / "results/annotation/bom_ss5_secretome.tsv"))
    ap.add_argument("--gtf", default=str(ROOT / "refs/BOM_ss5/BOM_ss5_genomic.gtf"))
    ap.add_argument("--outdir", default=str(ROOT / "results/wgcna"))
    ap.add_argument("--min-term", type=int, default=3)
    ap.add_argument("--fdr", type=float, default=0.05)
    args = ap.parse_args()

    import re
    gene_of = {}
    for line in open(args.gtf):
        if line.startswith("#"): continue
        p = line.split("\t")
        if len(p) < 9 or p[2] != "CDS": continue
        g = re.search(r'gene_id "([^"]+)"', p[8]); pr = re.search(r'protein_id "([^"]+)"', p[8])
        if g and pr: gene_of.setdefault(pr.group(1), g.group(1))

    cazy = {}
    for r in csv.DictReader(open(args.cazy), delimiter="\t"):
        g = gene_of.get(r["protein_id"])
        if g: cazy[g] = r["cazy_class"]
    secreted = set()
    for r in csv.DictReader(open(args.secretome), delimiter="\t"):
        if r["secreted"] == "1":
            g = gene_of.get(r["protein_id"])
            if g: secreted.add(g)

    modules = collections.defaultdict(list)
    for r in csv.DictReader(open(args.modules)):
        modules[r["module"]].append(r["gene"])

    g2t, t2g, tname = load_terms(args.terms)
    network = {g for gs in modules.values() for g in gs}
    background = {g for g in network if g in g2t}
    N = len(background)
    print(f"network genes: {len(network):,}; with annotation (background): {N:,} "
          f"({100*N/len(network):.0f}%)")

    # tissue association from the FDR-corrected module-trait table
    assoc = {}
    mtc = ROOT / "results/wgcna/module_trait_correlation.csv"
    mtp = ROOT / "results/wgcna/module_trait_pvalue.csv"
    if mtc.exists():
        cor = list(csv.reader(open(mtc))); pv = list(csv.reader(open(mtp)))
        traits = cor[0][1:]
        flat = []
        for cr, pr in zip(cor[1:], pv[1:]):
            for j, t in enumerate(traits):
                flat.append((cr[0].replace("ME", ""), t, float(cr[j+1]), float(pr[j+1])))
        rej, q, *_ = multipletests([x[3] for x in flat], method="fdr_bh")
        for (m, t, r, _p), sig, qq in zip(flat, rej, q):
            if sig and abs(r) > 0.7:
                assoc.setdefault(m, []).append((t, r, qq))

    rows, titles = [], []
    for mod, genes in sorted(modules.items(), key=lambda kv: -len(kv[1])):
        if mod == "grey":
            continue
        mg = {g for g in genes if g in background}
        if len(mg) < 5:
            continue
        best = {}
        for onto in ONTOLOGIES:
            cand = []
            for (o, tid), gs in t2g.items():
                if o != onto: continue
                gb = gs & background
                k = len(gs & mg)
                if k < 2 or len(gb) < args.min_term: continue
                p = hypergeom.sf(k - 1, N, len(gb), len(mg))
                cand.append((p, tid, tname[(o, tid)], k, len(gb)))
            if not cand: continue
            cand.sort()
            rej, q, *_ = multipletests([c[0] for c in cand], method="fdr_bh")
            sig = [(c, qq) for c, qq, rj in zip(cand, q, rej) if rj and qq < args.fdr]
            for c, qq in sig[:15]:
                rows.append({"module": mod, "module_size": len(genes), "annotated": len(mg),
                             "ontology": onto, "term_id": c[1], "term_name": c[2],
                             "genes_in_module": c[3], "genes_in_background": c[4],
                             "p": f"{c[0]:.3g}", "fdr": f"{qq:.3g}"})
            if sig: best[onto] = [s[0] for s in sig]

        # ---- readable title ----
        def pick(onto, n=2):
            out = []
            for c in best.get(onto, []):
                nm = c[2]
                if any(v.lower() == nm.lower() or v.lower() in nm.lower() for v in VAGUE) and len(out) == 0 and len(best.get(onto, [])) > 1:
                    continue
                out.append(nm)
                if len(out) == n: break
            return out

        # GO_CC is a legitimate last resort: several modules enrich only for localisation
        # (e.g. red -> fungal-type vacuole / translation initiation complexes). Omitting it
        # made those look unenriched and pushed them to the CAZyme fallback.
        bits = pick("GO_BP", 2) or pick("KEGG", 2) or pick("GO_MF", 2) or pick("GO_CC", 2)
        loc = pick("GO_CC", 1)
        n_caz = sum(1 for g in genes if g in cazy)
        n_sec = sum(1 for g in genes if g in secreted)
        tis = assoc.get(mod, [])
        tis_s = "; ".join(f"{t} r={r:+.2f}" for t, r, _q in tis)

        # A CAZyme count is a composition descriptor, not evidence of enrichment. Modules
        # with nothing significant are labelled as such rather than given a functional name
        # they have not earned -- turquoise (n=490) and blue (n=401) are genuinely unenriched,
        # most likely because only 59% of network genes carry any annotation at all.
        if bits:
            title = " / ".join(bits[:2])
        else:
            title = "no enriched terms (unannotated-dominated)"
        extra = []
        if loc: extra.append(loc[0])
        if n_sec: extra.append(f"{n_sec} secreted")
        if n_caz: extra.append(f"{n_caz} CAZy")
        titles.append({
            "module": mod, "size": len(genes), "annotated": len(mg),
            "title": title, "compartment": loc[0] if loc else "",
            "n_secreted": n_sec, "n_cazyme": n_caz,
            "tissue_association_FDR": tis_s,
            "readable_name": f"{mod}: {title}" + (f" [{', '.join(extra)}]" if extra else "")
                             + (f" -- {tis_s}" if tis_s else ""),
        })

    out = pathlib.Path(args.outdir)
    with (out / "module_enrichment.tsv").open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t"); w.writeheader(); w.writerows(rows)
    with (out / "module_titles.tsv").open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(titles[0]), delimiter="\t"); w.writeheader(); w.writerows(titles)

    n_enr = len({r["module"] for r in rows})
    print(f"modules tested: {len(titles)}; with >=1 enriched term at FDR<{args.fdr}: {n_enr}")
    print(f"\n=== modules with a tissue association surviving FDR ===")
    for t in titles:
        if t["tissue_association_FDR"]:
            print(f"  {t['readable_name']}")
    print(f"\nwrote {out/'module_enrichment.tsv'} and {out/'module_titles.tsv'}")


if __name__ == "__main__":
    sys.exit(main())
