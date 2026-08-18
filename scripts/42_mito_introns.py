#!/usr/bin/env python3
"""P4 -- catalogue mitochondrial introns from split-read evidence and measure splicing.

Junctions are taken from the N operations of the HISAT2 CIGAR strings on the mitochondrial
contig, pooled across all 16 libraries. A junction is retained only if it is seen in at least
three libraries with at least ten supporting reads in total, which removes the long tail of
one-off misalignments (each library individually shows 51-74 raw junctions, most singletons).

Retained junctions are then intersected with group I/II intron predictions from Rfam
covariance-model search, so the split-read evidence and the structural prediction are
independent lines of support for the same features.

Splicing efficiency per junction = spliced reads / (spliced + reads crossing the donor site
without a gap). It is reported per library with the full replicate spread, and deliberately
NOT contrasted between tissues: the within-tissue range spans 20-34 fold, which is wider than
any between-tissue difference these libraries could resolve.
"""
import argparse, collections, csv, pathlib, re, subprocess, sys

ROOT = pathlib.Path(__file__).resolve().parents[1]
CONTIG = "CM148777.1"
SAM = str(ROOT / "envs/root/envs/env-hisat/bin/samtools")


def junctions_and_depth(bam, contig):
    """Return {(start,end): spliced_reads} and {pos: reads crossing without a gap}."""
    juncs = collections.Counter()
    unspliced = collections.Counter()
    p = subprocess.run([SAM, "view", bam, contig], capture_output=True, text=True)
    for line in p.stdout.splitlines():
        f = line.split("\t")
        if len(f) < 6:
            continue
        pos, cig = int(f[3]), f[5]
        ref = pos
        blocks = []
        for num, op in re.findall(r"(\d+)([MIDNSHP=X])", cig):
            n = int(num)
            if op in "MD=X":
                blocks.append((ref, ref + n)); ref += n
            elif op == "N":
                juncs[(ref, ref + n)] += 1; ref += n
        if "N" not in cig:
            for a, b in blocks:
                unspliced[(a, b)] += 1
    return juncs, unspliced


def spans(unspliced, site, flank=8):
    """Reads whose aligned block spans a site with flanking bases on both sides."""
    return sum(c for (a, b), c in unspliced.items() if a + flank <= site <= b - flank)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bamdir", default=str(ROOT / "bam_BOM_ss5"))
    ap.add_argument("--cm", default=str(ROOT / "results/mito/introns_cm.tbl"))
    ap.add_argument("--outdir", default=str(ROOT / "results/mito"))
    ap.add_argument("--min-libs", type=int, default=3)
    ap.add_argument("--min-reads", type=int, default=10)
    args = ap.parse_args()
    out = pathlib.Path(args.outdir); out.mkdir(parents=True, exist_ok=True)

    tis = {x["sample_name"]: x["Factor Value[Tissue]"]
           for x in csv.DictReader((ROOT / "metadata/runsheet.csv").open())}

    per_lib, per_lib_unspl, raw_counts = {}, {}, {}
    for bam in sorted(pathlib.Path(args.bamdir).glob("*.primary.bam")):
        s = bam.name.replace(".primary.bam", "")
        j, u = junctions_and_depth(str(bam), CONTIG)
        per_lib[s] = j; per_lib_unspl[s] = u; raw_counts[s] = len(j)
        print(f"  {s.split('_sample_')[-1]:<4} raw junctions {len(j):>4}  spliced reads {sum(j.values()):>7,}", flush=True)

    allj = collections.Counter()
    nlib = collections.Counter()
    for s, j in per_lib.items():
        for k, v in j.items():
            allj[k] += v; nlib[k] += 1
    keep = [k for k in allj if nlib[k] >= args.min_libs and allj[k] >= args.min_reads]
    keep.sort(key=lambda k: -allj[k])
    print(f"\nraw junctions pooled: {len(allj):,}")
    print(f"retained (>={args.min_libs} libraries, >={args.min_reads} reads): {len(keep)}")

    # Rfam covariance-model intron predictions
    cm = []
    if pathlib.Path(args.cm).exists():
        for line in open(args.cm):
            if line.startswith("#"):
                continue
            f = line.split()
            if len(f) > 15:
                a, b = int(f[7]), int(f[8])
                cm.append((min(a, b), max(a, b), f[2], float(f[15])))
    def cm_hit(s, e):
        for a, b, name, ev in cm:
            if not (e < a or s > b):
                return f"{name} (E={ev:.1e})"
        return ""

    rows = []
    for (s, e) in keep:
        length = e - s
        libs = {k: per_lib[k].get((s, e), 0) for k in per_lib}
        effs = []
        for k in per_lib:
            sp = libs[k]
            un = spans(per_lib_unspl[k], s)
            if sp + un >= 5:
                effs.append(sp / (sp + un))
        rows.append(dict(start=s, end=e, length=length, total_reads=allj[(s, e)],
                         libraries=nlib[(s, e)], rfam=cm_hit(s, e),
                         mean_efficiency=round(sum(effs)/len(effs), 3) if effs else "",
                         min_efficiency=round(min(effs), 3) if effs else "",
                         max_efficiency=round(max(effs), 3) if effs else "",
                         n_libs_scored=len(effs),
                         **{f"reads_{k.split('_sample_')[-1]}": libs[k] for k in sorted(per_lib)}))

    with (out / "intron_catalogue.csv").open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)

    print(f"\n{'start':>7} {'end':>7} {'len':>6} {'reads':>8} {'libs':>5}  {'eff (min-max)':<18} rfam")
    for r in rows[:14]:
        eff = (f"{r['mean_efficiency']} ({r['min_efficiency']}-{r['max_efficiency']})"
               if r["mean_efficiency"] != "" else "n/a")
        print(f"{r['start']:>7} {r['end']:>7} {r['length']:>6} {r['total_reads']:>8,} "
              f"{r['libraries']:>5}  {eff:<18} {r['rfam']}")
    n_cm = sum(1 for r in rows if r["rfam"])
    print(f"\njunctions overlapping an Rfam intron prediction: {n_cm}/{len(rows)}")
    print(f"wrote {out/'intron_catalogue.csv'}")


if __name__ == "__main__":
    sys.exit(main())
