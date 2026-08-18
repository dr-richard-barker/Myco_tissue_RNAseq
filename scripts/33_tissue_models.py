#!/usr/bin/env python3
"""Phase 6b -- extract context-specific metabolic models for each tissue with RIPTiDe.

RIPTiDe weights reactions by transcript abundance and finds a flux distribution that is both
feasible and maximally consistent with the transcriptome, then samples the resulting solution
space. It is well suited to a small number of conditions and does not require an arbitrary
expression threshold, which matters here because absolute expression is low.

Scope, stated plainly: contextualisation uses aggregate expression across the whole network
rather than individual marker calls, so it tolerates per-gene noise better than the marker
analysis does. Even so, the exudophore and nodule models rest on the two tissues whose
signatures survived replicate-level and FDR-corrected testing; the two mycelial models are
included as a baseline and should be read with that asymmetry in mind.

Model gene identifiers are protein accessions (KAN*) while counts are keyed on locus tags
(AAD021_*); the mapping comes from the CDS records of the stock GTF.
"""
import argparse, collections, csv, json, pathlib, re, sys, warnings
warnings.filterwarnings("ignore")
import cobra, riptide

ROOT = pathlib.Path(__file__).resolve().parents[1]


def gene_to_protein(gtf):
    m = {}
    for line in open(gtf):
        if line.startswith("#"):
            continue
        p = line.split("\t")
        if len(p) < 9 or p[2] != "CDS":
            continue
        g = re.search(r'gene_id "([^"]+)"', p[8]); pr = re.search(r'protein_id "([^"]+)"', p[8])
        if g and pr:
            m.setdefault(g.group(1), pr.group(1))
    return m


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--model", default=str(ROOT / "models/BOM_ss5_gapfilled.xml"))
    ap.add_argument("--counts", default=str(ROOT / "results/normalised_all_genes.csv"))
    ap.add_argument("--gtf", default=str(ROOT / "refs/BOM_ss5/BOM_ss5_genomic.gtf"))
    ap.add_argument("--outdir", default=str(ROOT / "models/tissue"))
    ap.add_argument("--samples", type=int, default=500, help="flux samples per tissue")
    args = ap.parse_args()

    out = pathlib.Path(args.outdir); out.mkdir(parents=True, exist_ok=True)
    model = cobra.io.read_sbml_model(args.model)
    print(f"base model: {len(model.reactions)} reactions, {len(model.genes)} genes", flush=True)

    g2p = gene_to_protein(args.gtf)
    model_genes = {g.id for g in model.genes}

    with open(args.counts) as fh:
        r = csv.reader(fh); hdr = next(r); samples = hdr[1:]
        counts = {row[0]: [float(x) for x in row[1:]] for row in r}

    tis = {x["sample_name"]: x["Factor Value[Tissue]"]
           for x in csv.DictReader((ROOT / "metadata/runsheet.csv").open())}
    groups = collections.defaultdict(list)
    for i, s in enumerate(samples):
        groups[tis[s]].append(i)

    summary = []
    for tissue in sorted(groups):
        idx = groups[tissue]
        transcriptome, mapped = {}, 0
        for gene, vals in counts.items():
            prot = g2p.get(gene)
            if prot and prot in model_genes:
                transcriptome[prot] = sum(vals[i] for i in idx) / len(idx)
                mapped += 1
        print(f"\n=== {tissue} (n={len(idx)}) : {mapped} genes mapped into the model ===", flush=True)
        if mapped < 20:
            print("  too few mapped genes; skipping", flush=True)
            continue
        try:
            res = riptide.contextualize(model=model, transcriptome=transcriptome,
                                        samples=args.samples, silent=True)
        except Exception as e:
            print(f"  RIPTiDe failed: {type(e).__name__}: {str(e)[:120]}", flush=True)
            continue

        cm = res.model
        tag = tissue.replace(" ", "_")
        cobra.io.write_sbml_model(cm, str(out / f"BOM_ss5_{tag}.xml"))
        flux = res.flux_samples
        flux.to_csv(out / f"flux_{tag}.csv", index=False)
        obj = float(res.fraction_of_optimum) if hasattr(res, "fraction_of_optimum") else float("nan")
        print(f"  reactions retained: {len(cm.reactions)} / {len(model.reactions)}", flush=True)
        print(f"  metabolites: {len(cm.metabolites)}   genes: {len(cm.genes)}", flush=True)
        summary.append(dict(tissue=tissue, n=len(idx), mapped_genes=mapped,
                            reactions=len(cm.reactions), metabolites=len(cm.metabolites),
                            genes=len(cm.genes), flux_samples=flux.shape[0]))

    if summary:
        with (out / "tissue_model_summary.csv").open("w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(summary[0])); w.writeheader(); w.writerows(summary)
        print(f"\nwrote {out}/tissue_model_summary.csv")
    return 0


if __name__ == "__main__":
    sys.exit(main())
