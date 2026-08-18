#!/usr/bin/env python3
"""Phase 5/6 -- turn DIAMOND hits against Swiss-Prot into a functional annotation table.

eggNOG-mapper is the usual tool here, but eggnogdb.embl.de is unreachable from this machine
(DNS failure; UniProt, EBI, NCBI and GenomeNet all resolve fine). Swiss-Prot is a reasonable
substitute for the two things downstream work actually needs: EC numbers to seed the
genome-scale metabolic reconstruction, and GO terms plus protein names for the per-tissue
transcriptional model.

Reads the .dat flatfile only for accessions that were actually hit, so the 700 MB archive is
streamed once and nothing is held in memory beyond the hit set.
"""
import argparse
import collections
import csv
import gzip
import pathlib
import sys

ROOT = pathlib.Path(__file__).resolve().parents[1]


def best_hits(tsv, min_pident, max_evalue):
    """One best hit per query, by bitscore."""
    best = {}
    with open(tsv) as fh:
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) < 6:
                continue
            q, s, pident, _len, evalue, bits = p[0], p[1], float(p[2]), p[3], float(p[4]), float(p[5])
            if pident < min_pident or evalue > max_evalue:
                continue
            if q not in best or bits > best[q][2]:
                # sseqid looks like sp|P12345|NAME_ORG
                acc = s.split("|")[1] if "|" in s else s
                best[q] = (acc, pident, bits, evalue)
    return best


def scan_dat(dat_gz, wanted):
    """Stream the Swiss-Prot flatfile, pulling EC / GO / name for wanted accessions."""
    info = {}
    acc = None
    keep = False
    rec = None
    with gzip.open(dat_gz, "rt", errors="replace") as fh:
        for line in fh:
            if line.startswith("AC "):
                if acc is None:
                    for a in line[5:].replace(" ", "").rstrip(";\n").split(";"):
                        if a in wanted:
                            acc, keep = a, True
                            rec = {"ec": set(), "go": set(), "name": "", "kegg": set()}
                            break
                    else:
                        acc, keep = "?", False
            elif line.startswith("DE ") and keep:
                if "EC=" in line:
                    for frag in line.split("EC=")[1:]:
                        rec["ec"].add(frag.split(";")[0].strip().rstrip(","))
                if "RecName: Full=" in line and not rec["name"]:
                    rec["name"] = line.split("RecName: Full=")[1].split(";")[0].split("{")[0].strip()
            elif line.startswith("DR ") and keep:
                if line.startswith("DR   GO;"):
                    rec["go"].add(line.split(";")[1].strip())
                elif line.startswith("DR   KEGG;"):
                    rec["kegg"].add(line.split(";")[1].strip())
            elif line.startswith("//"):
                if keep and acc and acc != "?":
                    info[acc] = rec
                acc, keep, rec = None, False, None
    return info


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--hits", default=str(ROOT / "results" / "annotation" / "pc915_vs_sprot.tsv"))
    ap.add_argument("--dat", default=str(ROOT / "refs" / "uniprot" / "uniprot_sprot.dat.gz"))
    ap.add_argument("--out", default=str(ROOT / "results" / "annotation" / "pc915_functional.tsv"))
    ap.add_argument("--min-pident", type=float, default=30.0)
    ap.add_argument("--max-evalue", type=float, default=1e-5)
    args = ap.parse_args()

    if not pathlib.Path(args.hits).exists():
        sys.exit(f"missing {args.hits} -- run diamond blastp first")

    best = best_hits(args.hits, args.min_pident, args.max_evalue)
    print(f"proteins with a Swiss-Prot hit: {len(best):,}")
    wanted = {v[0] for v in best.values()}
    print(f"distinct Swiss-Prot accessions to look up: {len(wanted):,}")

    info = scan_dat(args.dat, wanted)
    print(f"accessions resolved in flatfile: {len(info):,}")

    out = pathlib.Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    n_ec = 0
    with out.open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["protein_id", "sprot_acc", "pident", "evalue", "protein_name", "EC", "GO", "KEGG"])
        for q, (acc, pident, _bits, evalue) in sorted(best.items()):
            rec = info.get(acc, {})
            ecs = ";".join(sorted(rec.get("ec", ())))
            if ecs:
                n_ec += 1
            w.writerow([q, acc, f"{pident:.1f}", f"{evalue:.1e}", rec.get("name", ""),
                        ecs, ";".join(sorted(rec.get("go", ()))), ";".join(sorted(rec.get("kegg", ())))])

    print(f"proteins with at least one EC number: {n_ec:,}")
    ecset = collections.Counter()
    for q, (acc, *_rest) in best.items():
        for e in info.get(acc, {}).get("ec", ()):
            ecset[e] += 1
    print(f"distinct EC numbers: {len(ecset):,}  (these seed the draft GEM)")
    print(f"wrote {out}")


if __name__ == "__main__":
    sys.exit(main())
