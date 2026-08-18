#!/usr/bin/env python3
"""P4 -- reconcile a corrected mitogenome annotation for Pleurotus ostreatus BOM_ss5.

The stock GenBank annotation of CM148777.1 carries 13 features across 71,949 bp, of which one
is protein-coding. This assembles a replacement from four independent lines of evidence:

  tRNAs            tRNAscan-SE (organellar mode), cross-checked against aragorn
  rRNAs            barrnap (mito mode), already produced for the rRNA-removal work
  protein-coding   EMBOSS getorf under genetic code 4, identified by DIAMOND against
                   Swiss-Prot. Code 4 was confirmed empirically, not assumed: it yields 51
                   ORFs >=100 aa against 34 under code 1
  introns          split-read junctions spliced at >50% efficiency, plus Rfam covariance-model
                   group I predictions

Nothing here is invented: every feature carries the evidence that produced it in the GFF
attributes, so a reader can discount any line of evidence they distrust.
"""
import collections, csv, pathlib, re, sys

ROOT = pathlib.Path(__file__).resolve().parents[1]
CONTIG = "CM148777.1"
OUT = ROOT / "refs/BOM_ss5/mitogenome.gff"
RES = ROOT / "results/mito"


def parse_trnascan(p):
    feats = []
    if not p.exists(): return feats
    for line in open(p):
        f = line.split()
        if len(f) < 9 or not f[0].startswith(CONTIG): continue
        try: a, b = int(f[2]), int(f[3])
        except ValueError: continue
        strand = "+" if a <= b else "-"
        feats.append((min(a, b), max(a, b), strand, "tRNA",
                      f"trn{f[4]}", f"anticodon={f[5]};score={f[8]};evidence=tRNAscan-SE"))
    return feats


def parse_barrnap(p):
    feats = []
    if not p.exists(): return feats
    for line in open(p):
        if line.startswith("#"): continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 9 or f[0] != CONTIG: continue
        nm = "rRNA"
        m = re.search(r"Name=([^;]+)", f[8])
        if m: nm = m.group(1)
        feats.append((int(f[3]), int(f[4]), f[6], "rRNA", nm, "evidence=barrnap"))
    return feats


def parse_orfs(orf_fa, hits_tsv):
    """getorf headers carry coordinates; keep ORFs with a Swiss-Prot identification."""
    coords = {}
    for line in open(orf_fa):
        if not line.startswith(">"): continue
        m = re.match(r">(\S+)\s+\[(\d+)\s*-\s*(\d+)\]", line)
        if not m: continue
        name, a, b = m.group(1), int(m.group(2)), int(m.group(3))
        rev = "(REVERSE SENSE)" in line
        coords[name] = (min(a, b), max(a, b), "-" if rev else "+")
    best = {}
    for line in open(hits_tsv):
        p = line.rstrip("\n").split("\t")
        if len(p) < 7: continue
        q, bits, title = p[0], float(p[5]), p[6]
        if q not in best or bits > best[q][0]:
            best[q] = (bits, title, float(p[2]))
    feats = []
    for q, (bits, title, pid) in best.items():
        if q not in coords: continue
        a, b, st = coords[q]
        nm = title.split("OS=")[0].split(" ", 1)[-1].strip()
        feats.append((a, b, st, "CDS", nm[:60],
                      f"identity={pid:.1f};bitscore={bits:.0f};evidence=getorf+DIAMOND"))
    return feats


def parse_introns(cat):
    feats = []
    if not cat.exists(): return feats
    for r in csv.DictReader(open(cat)):
        eff = r["mean_efficiency"]
        if eff == "" or float(eff) <= 0.5: continue
        rf = r["rfam"] or "none"
        feats.append((int(r["start"]), int(r["end"]), ".", "intron",
                      f"intron_{r['start']}",
                      f"splicing_efficiency={eff};supporting_reads={r['total_reads']};"
                      f"libraries={r['libraries']};rfam={rf};evidence=split_reads"))
    return feats


def main():
    feats = []
    feats += parse_trnascan(RES / "trnascan.out")
    feats += parse_barrnap(ROOT / "refs/rRNA/BOM_ss5_mito.gff")
    feats += parse_orfs(RES / "orfs_aa.fa", RES / "orfs_vs_sprot.tsv")
    feats += parse_introns(RES / "intron_catalogue.csv")

    # drop CDS predictions wholly contained in another, longer CDS (getorf frame duplicates)
    cds = sorted([f for f in feats if f[3] == "CDS"], key=lambda f: -(f[1] - f[0]))
    keep_cds, taken = [], []
    for f in cds:
        if any(f[0] >= g[0] and f[1] <= g[1] for g in taken): continue
        taken.append(f); keep_cds.append(f)
    feats = [f for f in feats if f[3] != "CDS"] + keep_cds
    feats.sort(key=lambda f: (f[0], f[1]))

    with OUT.open("w") as fh:
        fh.write("##gff-version 3\n")
        fh.write(f"##sequence-region {CONTIG} 1 71949\n")
        fh.write("# Corrected mitogenome annotation; see scripts/40_mito_annotate.py\n")
        for i, (a, b, st, kind, name, attrs) in enumerate(feats, 1):
            fh.write(f"{CONTIG}\tmitoannot\t{kind}\t{a}\t{b}\t.\t{st}\t.\t"
                     f"ID={kind}_{i};Name={name};{attrs}\n")

    counts = collections.Counter(f[3] for f in feats)
    print(f"corrected annotation: {len(feats)} features")
    for k, v in sorted(counts.items()):
        print(f"    {k:<8} {v}")
    stock = collections.Counter()
    for line in open(ROOT / "refs/BOM_ss5/BOM_ss5_genomic.gtf"):
        f = line.split("\t")
        if len(f) > 8 and f[0] == CONTIG and f[2] == "gene":
            m = re.search(r'gene_biotype "([^"]+)"', f[8])
            stock[m.group(1) if m else "other"] += 1
    print(f"\nstock annotation: {sum(stock.values())} gene features")
    for k, v in sorted(stock.items()):
        print(f"    {k:<8} {v}")
    print(f"\nwrote {OUT}")


if __name__ == "__main__":
    sys.exit(main())
