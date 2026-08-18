#!/usr/bin/env python3
"""Phase 2b -- extend gene 3' ends so 3'-tag reads land inside a counted feature.

These libraries are 3'-end tag counting: reads pile up at the polyadenylation site, which
sits downstream of the stop codon. Even PC9.15, the best-annotated candidate, has a median
annotated 3'UTR of only 78 bp and leaves ~25% of transcripts with none at all. The
consequence is measurable: featureCounts assigns only ~15% of rDNA-depleted reads to genes
while 30-78% land in "no feature".

Extension rules:
  * extend each gene's 3' end (strand-aware) by up to --max bp;
  * cap at the distance to the next annotated feature on EITHER strand minus --gap, so a
    compact fungal genome does not get read-through cross-assignment between neighbours;
  * never extend past the end of the contig.

Emits a GTF whose `exon` records carry the extension, since featureCounts counts on -t exon.
Only the last exon of each transcript is lengthened, so internal structure is untouched.
"""
import argparse
import collections
import gzip
import pathlib
import sys


def openmaybe(p):
    return gzip.open(p, "rt") if str(p).endswith(".gz") else open(p)


def attr(field, key):
    marker = key + ' "'
    i = field.find(marker)
    if i < 0:
        return None
    j = field.find('"', i + len(marker))
    return field[i + len(marker):j]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("gtf")
    ap.add_argument("out")
    ap.add_argument("--max", type=int, default=500, help="max 3' extension (bp)")
    ap.add_argument("--gap", type=int, default=50, help="buffer to keep before next feature")
    args = ap.parse_args()

    lines = []
    # gene span and strand, plus every feature interval per contig for the neighbour cap
    gene_span = {}
    occupied = collections.defaultdict(list)
    contig_len = collections.defaultdict(int)

    with openmaybe(args.gtf) as fh:
        for raw in fh:
            if raw.startswith("#"):
                continue
            p = raw.rstrip("\n").split("\t")
            if len(p) < 9:
                continue
            lines.append(p)
            start, end = int(p[3]), int(p[4])
            contig_len[p[0]] = max(contig_len[p[0]], end)
            if p[2] == "gene":
                gid = attr(p[8], "gene_id")
                if gid:
                    gene_span[gid] = [p[0], start, end, p[6]]
            if p[2] in ("gene", "exon"):
                occupied[p[0]].append((start, end))

    for c in occupied:
        occupied[c].sort()

    def next_boundary(contig, pos, strand, gid_span):
        """Nearest feature edge downstream of pos, ignoring the gene's own span."""
        gs, ge = gid_span
        best = None
        for s, e in occupied[contig]:
            if s >= gs and e <= ge:
                continue  # the gene's own features
            if strand == "+":
                if s > pos and (best is None or s < best):
                    best = s
            else:
                if e < pos and (best is None or e > best):
                    best = e
        return best

    # decide the extension per gene
    ext = {}
    for gid, (contig, gs, ge, strand) in gene_span.items():
        if strand == "+":
            limit = args.max
            nb = next_boundary(contig, ge, "+", (gs, ge))
            if nb is not None:
                limit = min(limit, max(0, nb - args.gap - ge))
            limit = min(limit, max(0, contig_len[contig] - ge))
        else:
            limit = args.max
            nb = next_boundary(contig, gs, "-", (gs, ge))
            if nb is not None:
                limit = min(limit, max(0, gs - args.gap - nb))
            limit = min(limit, max(0, gs - 1))
        ext[gid] = limit

    # find the terminal exon of each transcript so only that one grows
    tx_bounds = {}
    for p in lines:
        if p[2] != "exon":
            continue
        tid = attr(p[8], "transcript_id")
        if not tid:
            continue
        s, e = int(p[3]), int(p[4])
        b = tx_bounds.setdefault(tid, [s, e])
        b[0] = min(b[0], s)
        b[1] = max(b[1], e)

    n_ext = 0
    total_bp = 0
    out = pathlib.Path(args.out)
    with out.open("w") as fh:
        for p in lines:
            gid = attr(p[8], "gene_id")
            e = ext.get(gid, 0)
            if e > 0 and p[2] in ("gene", "transcript", "exon"):
                strand = p[6]
                grow = False
                if p[2] in ("gene", "transcript"):
                    grow = True
                else:
                    tid = attr(p[8], "transcript_id")
                    b = tx_bounds.get(tid)
                    if b:
                        grow = (int(p[4]) == b[1]) if strand == "+" else (int(p[3]) == b[0])
                if grow:
                    if strand == "+":
                        p[4] = str(int(p[4]) + e)
                    else:
                        p[3] = str(max(1, int(p[3]) - e))
                    if p[2] == "gene":
                        n_ext += 1
                        total_bp += e
            fh.write("\t".join(p) + "\n")

    nonzero = [v for v in ext.values() if v > 0]
    capped = sum(1 for v in ext.values() if 0 < v < args.max)
    print(f"genes: {len(ext):,}")
    print(f"extended: {n_ext:,} ({100 * n_ext / max(len(ext), 1):.1f}%)")
    print(f"  at full {args.max} bp : {sum(1 for v in ext.values() if v == args.max):,}")
    print(f"  capped by neighbour  : {capped:,}")
    print(f"  no room (0 bp)       : {sum(1 for v in ext.values() if v == 0):,}")
    if nonzero:
        nonzero.sort()
        print(f"  median extension     : {nonzero[len(nonzero) // 2]} bp")
    print(f"total added: {total_bp / 1e3:.1f} kb")
    print(f"wrote {out}")


if __name__ == "__main__":
    sys.exit(main())
