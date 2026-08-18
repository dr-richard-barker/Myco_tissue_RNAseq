#!/usr/bin/env python3
"""Phase 2 -- build the final counting annotation for a chosen reference.

Combines three things the stock GTFs lack:

1. 3' extension (from 07_extend_3p.py) so 3'-tag reads fall inside a counted feature.

2. Nuclear rRNA features from barrnap. None of the candidate GTFs annotate rDNA, so rRNA
   reads land in "no feature" and are invisible to the GL-DPPD-7101-G rRNA-removed track,
   which works by dropping rRNA-biotype features from the count matrix.

3. Coverage-derived rRNA blocks, supplied as a TSV (refs/rRNA_regions_<label>.tsv). barrnap
   calls only fragments of the real transcribed regions; the blocks come from clustering
   pooled intergenic reads and are confirmed by the barrnap calls falling inside them.
   For PC9.15 this is the mitochondrial rRNA (60.9% of pooled aligned reads); for BOM_ss5 it
   is that plus the resolved tandem nuclear rDNA array (~71% combined).

Everything added carries gene_biotype "rRNA" so the Phase 4 rRNA-removed track drops it by
biotype, exactly as GeneLab does. Blocks are emitted on a single strand: annotating both
makes every overlapping read Unassigned_Ambiguity under unstranded counting, which silently
discarded 27k-102k reads per sample when this was first written.
"""
import argparse
import pathlib
import sys


def gtf_record(contig, start, end, strand, feature, gene_id, name):
    attrs = (f'gene_id "{gene_id}"; transcript_id "{gene_id}.1"; '
             f'gene_name "{name}"; gene_biotype "rRNA";')
    return f"{contig}\tcustom\t{feature}\t{start}\t{end}\t.\t{strand}\t.\t{attrs}"


def load_regions(path):
    out = []
    if not path or not pathlib.Path(path).exists():
        return out
    for raw in open(path):
        if raw.startswith("#") or not raw.strip():
            continue
        p = raw.rstrip("\n").split("\t")
        if len(p) < 4:
            continue
        out.append((p[0], int(p[1]), int(p[2]), p[3]))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("base_gtf", help="3'-extended GTF from 07_extend_3p.py")
    ap.add_argument("barrnap_bed", help="rRNA BED from 04_rrna_ref.sh")
    ap.add_argument("out")
    ap.add_argument("--regions", help="TSV of coverage-derived rRNA blocks")
    args = ap.parse_args()

    blocks = load_regions(args.regions)
    block_contigs = {b[0] for b in blocks}

    out = pathlib.Path(args.out)
    added_barrnap, added_blocks = 0, 0
    with out.open("w") as fh:
        for line in open(args.base_gtf):
            fh.write(line)

        # barrnap loci, skipping contigs already covered by an explicit block (the blocks
        # are broader and would otherwise double-annotate the same reads)
        for raw in open(args.barrnap_bed):
            p = raw.rstrip("\n").split("\t")
            if len(p) < 6:
                continue
            contig, start, end, name, _, strand = p[:6]
            if contig in block_contigs:
                continue
            gid = f"rRNA_{contig}_{start}"
            for feat in ("gene", "transcript", "exon"):
                fh.write(gtf_record(contig, int(start) + 1, int(end), strand, feat, gid, name) + "\n")
            added_barrnap += 1

        for contig, start, end, name in blocks:
            gid = f"rRNA_{contig}_{start}"
            for feat in ("gene", "transcript", "exon"):
                fh.write(gtf_record(contig, start, end, "+", feat, gid, name) + "\n")
            added_blocks += 1

    print(f"barrnap loci added        : {added_barrnap}")
    print(f"coverage-derived blocks   : {added_blocks}")
    print(f"wrote {out}")


if __name__ == "__main__":
    sys.exit(main())
