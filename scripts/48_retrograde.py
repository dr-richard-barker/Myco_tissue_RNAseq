#!/usr/bin/env python3
"""P5 -- predicted nuclear-encoded mitochondrial proteome, and mito-nuclear coordination.

MitoCarta3.0 is a mammalian inventory, so orthology transfer to a basidiomycete is imperfect
in both directions: fungal mitochondrial proteins with no mammalian counterpart (alternative
oxidase being the standard example) are invisible here. Reciprocal best hits are used rather
than one-way best hits specifically to keep the predicted set conservative; a permissive screen
would return implausibly many.

Swiss-Prot's own SUBCELLULAR LOCATION annotations are screened in parallel as a fungal-aware
complement, since Swiss-Prot contains curated yeast and filamentous-fungal entries that
MitoCarta by construction does not.

The retrograde section reports whether mtDNA-encoded and nuclear-encoded OXPHOS expression
covary across tissues. Covariation is consistent with coordination; it does not demonstrate
signalling, and is not presented as doing so.
"""
import collections, csv, gzip, os, pathlib, re, subprocess, sys

ROOT = pathlib.Path(__file__).resolve().parents[1]
MM = str(ROOT / "envs/bin/micromamba")
ENV = dict(os.environ, MAMBA_ROOT_PREFIX=str((ROOT / "envs/root").resolve()))
OUT = ROOT / "results/mitocarta"; OUT.mkdir(parents=True, exist_ok=True)


def dia(q, d, o, extra=()):
    subprocess.run([MM, "run", "-n", "env-annot", "diamond", "blastp", "-q", str(q), "-d", str(d),
                    "-o", str(o), "--threads", "8", "--evalue", "1e-5", "--max-target-seqs", "1",
                    "--quiet", "--outfmt", "6", "qseqid", "sseqid", "pident", "evalue", "bitscore",
                    *extra], env=ENV, capture_output=True)


def best(path):
    b = {}
    for line in open(path):
        f = line.rstrip("\n").split("\t")
        if len(f) < 5: continue
        if f[0] not in b or float(f[4]) > b[f[0]][1]:
            b[f[0]] = (f[1], float(f[4]), float(f[2]))
    return b


def main():
    nuc = ROOT / "refs/BOM_ss5/BOM_ss5_protein.faa"
    mc = ROOT / "refs/mitocarta/Human.MitoCarta3.0.fasta"

    fwd = OUT / "nuc_vs_mitocarta.tsv"
    rev = OUT / "mitocarta_vs_nuc.tsv"
    if not fwd.exists():
        dia(nuc, ROOT / "refs/mitocarta/mitocarta", fwd)
    if not rev.exists():
        subprocess.run([MM, "run", "-n", "env-annot", "diamond", "makedb", "--in", str(nuc),
                        "-d", "/tmp/nucdb", "--quiet"], env=ENV, capture_output=True)
        dia(mc, "/tmp/nucdb", rev)

    F, R = best(fwd), best(rev)
    rbh = {q: v for q, v in F.items() if R.get(v[0], ("",))[0] == q}
    print(f"nuclear proteins with any MitoCarta hit : {len(F):,}")
    print(f"reciprocal best hits (predicted mito)   : {len(rbh):,}")

    # fungal-aware complement: Swiss-Prot curated mitochondrial localisation
    want = set()
    func = ROOT / "results/annotation/bom_ss5_functional.tsv"
    prot2acc = {}
    for r in csv.DictReader(func.open(), delimiter="\t"):
        if r.get("sprot_acc"):
            prot2acc[r["protein_id"]] = r["sprot_acc"]; want.add(r["sprot_acc"])
    mito_acc = set()
    acc = None; keep = False; is_mito = False
    with gzip.open(ROOT / "refs/uniprot/uniprot_sprot.dat.gz", "rt", errors="replace") as fh:
        for line in fh:
            if line.startswith("AC "):
                if acc is None:
                    for a in line[5:].replace(" ", "").rstrip(";\n").split(";"):
                        if a in want: acc, keep = a, True; break
                    else: acc, keep = "?", False
            elif keep and line.startswith("CC       ") and "Mitochond" in line:
                is_mito = True
            elif keep and line.startswith("CC   -!- SUBCELLULAR LOCATION") and "Mitochond" in line:
                is_mito = True
            elif line.startswith("//"):
                if keep and is_mito and acc != "?": mito_acc.add(acc)
                acc, keep, is_mito = None, False, False
    sp_mito = {p for p, a in prot2acc.items() if a in mito_acc}
    print(f"Swiss-Prot curated mitochondrial        : {len(sp_mito):,}")
    union = set(rbh) | sp_mito
    print(f"union (predicted mito-proteome)         : {len(union):,} "
          f"({100*len(union)/12521:.1f}% of the proteome)")

    with (OUT / "nuclear_mito_proteome.csv").open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["protein_id", "mitocarta_rbh", "mitocarta_target", "identity", "swissprot_mito"])
        for p in sorted(union):
            v = rbh.get(p)
            w.writerow([p, int(p in rbh), v[0] if v else "", f"{v[2]:.1f}" if v else "",
                        int(p in sp_mito)])
    print(f"wrote {OUT/'nuclear_mito_proteome.csv'}")
    return union


if __name__ == "__main__":
    main()


def retrograde():
    """Do mtDNA-encoded and nuclear-encoded mitochondrial genes covary across tissues?"""
    import statistics
    tis = {x["sample_name"]: x["Factor Value[Tissue]"]
           for x in csv.DictReader((ROOT / "metadata/runsheet.csv").open())}
    with (ROOT / "results/normalised_all_genes.csv").open() as fh:
        r = csv.reader(fh); hdr = next(r); samples = hdr[1:]
        expr = {row[0]: [float(x) for x in row[1:]] for row in r}
    grp = collections.defaultdict(list)
    for i, s in enumerate(samples): grp[tis[s]].append(i)
    tissues = sorted(grp)

    # nuclear mito proteome -> locus tags
    p2g = {}
    for line in (ROOT / "refs/BOM_ss5/BOM_ss5_genomic.gtf").open():
        if line.startswith("#"): continue
        f = line.split("\t")
        if len(f) < 9 or f[2] != "CDS": continue
        g = re.search(r'gene_id "([^"]+)"', f[8]); pr = re.search(r'protein_id "([^"]+)"', f[8])
        if g and pr: p2g.setdefault(pr.group(1), g.group(1))
    nuc = [p2g[r["protein_id"]] for r in csv.DictReader((OUT / "nuclear_mito_proteome.csv").open())
           if p2g.get(r["protein_id"]) in expr]
    print(f"\nnuclear mito genes with expression: {len(nuc):,}")

    # mtDNA-encoded expression, per tissue, from the P4 count matrix
    mt = {}
    with (ROOT / "results/mito/mito_counts.txt").open() as fh:
        for line in fh:
            if line.startswith("#"): continue
            f = line.rstrip("\n").split("\t")
            if f[0] == "Geneid":
                mtsamp = [pathlib.Path(x).name.split(".")[0] for x in f[6:]]; continue
            mt[f[0]] = [int(x) for x in f[6:]]
    mtgrp = collections.defaultdict(list)
    for i, s in enumerate(mtsamp): mtgrp[tis[s]].append(i)
    # restrict to CDS features
    cds = set()
    for line in (ROOT / "refs/BOM_ss5/mitogenome.gff").open():
        if line.startswith("#"): continue
        f = line.split("\t")
        if len(f) > 8 and f[2] == "CDS":
            cds.add(re.search(r"ID=([^;]+)", f[8]).group(1))
    mt_tot = {t: sum(sum(v[i] for i in mtgrp[t]) / len(mtgrp[t])
                     for k, v in mt.items() if k in cds) for t in tissues}
    nuc_tot = {t: sum(statistics.mean(expr[g][i] for i in grp[t]) for g in nuc) for t in tissues}

    print(f"\n{'tissue':<20}{'mtDNA-encoded':>16}{'nuclear mito':>16}{'ratio':>9}")
    for t in tissues:
        print(f"{t:<20}{mt_tot[t]:>16,.0f}{nuc_tot[t]:>16,.0f}{mt_tot[t]/max(nuc_tot[t],1):>9.3f}")

    x = [mt_tot[t] for t in tissues]; y = [nuc_tot[t] for t in tissues]
    mx, my = statistics.mean(x), statistics.mean(y)
    num = sum((a-mx)*(b-my) for a, b in zip(x, y))
    den = (sum((a-mx)**2 for a in x) * sum((b-my)**2 for b in y)) ** .5
    print(f"\ncorr(mtDNA-encoded, nuclear mito) across {len(tissues)} tissues = "
          f"{num/den if den else float('nan'):+.3f}")
    print("  n=4 tissues; this is descriptive, not a test.")

    with (OUT / "retrograde_summary.csv").open("w", newline="") as fh:
        w = csv.writer(fh); w.writerow(["tissue", "mtDNA_encoded", "nuclear_mito", "ratio"])
        for t in tissues:
            w.writerow([t, round(mt_tot[t], 1), round(nuc_tot[t], 1),
                        round(mt_tot[t]/max(nuc_tot[t], 1), 4)])
    print(f"wrote {OUT/'retrograde_summary.csv'}")


retrograde()
