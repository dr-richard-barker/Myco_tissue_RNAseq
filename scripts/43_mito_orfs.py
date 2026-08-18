#!/usr/bin/env python3
"""P4 -- free-standing ORFs of unknown function in the mitogenome, and whether they are expressed.

The Bipolaris sorokiniana mitogenome carries 52 free-standing ORFs of unknown function. This
asks the same question of P. ostreatus: which predicted ORFs have no characterised homologue,
and of those, which are actually transcribed.

An ORF is called "uncharacterised" only if it has no DIAMOND hit to Swiss-Prot at all. That is
a weaker statement than it sounds -- Swiss-Prot covers 41% of this proteome overall -- so the
set is a candidate list for follow-up, not a claim of novelty.
"""
import csv, pathlib, re, subprocess, sys, collections

ROOT = pathlib.Path(__file__).resolve().parents[1]
RES = ROOT / "results/mito"
SAM = str(ROOT / "envs/root/envs/env-hisat/bin/samtools")
CONTIG = "CM148777.1"


def orf_coords(fa):
    d = {}
    for line in open(fa):
        if not line.startswith(">"): continue
        m = re.match(r">(\S+)\s+\[(\d+)\s*-\s*(\d+)\]", line)
        if not m: continue
        a, b = int(m.group(2)), int(m.group(3))
        d[m.group(1)] = (min(a, b), max(a, b), "-" if "(REVERSE SENSE)" in line else "+")
    return d


def main():
    coords = orf_coords(RES / "orfs_aa.fa")
    hit = set()
    for line in open(RES / "orfs_vs_sprot.tsv"):
        hit.add(line.split("\t")[0])
    # annotated features to subtract (tRNA / rRNA / identified CDS)
    ann = []
    for line in open(ROOT / "refs/BOM_ss5/mitogenome.gff"):
        if line.startswith("#"): continue
        f = line.split("\t")
        if len(f) < 9 or f[2] == "intron": continue
        ann.append((int(f[3]), int(f[4])))

    def overlaps_annotated(a, b, frac=0.5):
        for x, y in ann:
            ov = min(b, y) - max(a, x)
            if ov > 0 and ov >= frac * (b - a):
                return True
        return False

    cand = {k: v for k, v in coords.items() if k not in hit and not overlaps_annotated(*v[:2])}
    print(f"ORFs predicted (>=100 aa, code 4): {len(coords)}")
    print(f"  with a Swiss-Prot hit          : {len(hit & set(coords))}")
    print(f"  no hit and not overlapping an annotated feature: {len(cand)}")

    # expression: reads overlapping each candidate, per library
    bams = sorted((ROOT / "bam_BOM_ss5").glob("*.primary.bam"))
    tis = {x["sample_name"]: x["Factor Value[Tissue]"]
           for x in csv.DictReader((ROOT / "metadata/runsheet.csv").open())}
    rows = []
    for name, (a, b, st) in sorted(cand.items(), key=lambda kv: kv[1][0]):
        per = {}
        for bam in bams:
            s = bam.name.replace(".primary.bam", "")
            r = subprocess.run([SAM, "view", "-c", str(bam), f"{CONTIG}:{a}-{b}"],
                               capture_output=True, text=True)
            per[s] = int(r.stdout.strip() or 0)
        tot = sum(per.values())
        det = sum(1 for v in per.values() if v >= 5)
        rows.append(dict(orf=name, start=a, end=b, strand=st, length_nt=b - a + 1,
                         total_reads=tot, libraries_detected=det,
                         **{f"reads_{k.split('_sample_')[-1]}": v for k, v in sorted(per.items())}))

    rows.sort(key=lambda r: -r["total_reads"])
    with (RES / "candidate_orfs.csv").open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)

    expressed = [r for r in rows if r["libraries_detected"] >= 8]
    print(f"\nexpressed in >=8 of 16 libraries: {len(expressed)}")
    print(f"{'ORF':<22}{'start':>8}{'end':>8}{'nt':>7}{'reads':>10}{'libs':>6}")
    for r in rows[:12]:
        print(f"{r['orf']:<22}{r['start']:>8}{r['end']:>8}{r['length_nt']:>7}"
              f"{r['total_reads']:>10,}{r['libraries_detected']:>6}")
    print(f"\nwrote {RES/'candidate_orfs.csv'}")


if __name__ == "__main__":
    sys.exit(main())
