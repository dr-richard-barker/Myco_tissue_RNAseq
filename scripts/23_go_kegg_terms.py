#!/usr/bin/env python3
"""Phase 5 -- build GO (CC/MF/BP) and KEGG pathway term tables for the annotated proteome.

The earlier annotation (13_annotate.py) kept only bare GO IDs. Enrichment needs the ontology
ASPECT and the human-readable term name, so Swiss-Prot is re-parsed here: its GO cross-refs
carry both, e.g.
    DR   GO; GO:0005576; C:extracellular region; IEA:UniProtKB-SubCell.
where the leading letter is C (cellular component), F (molecular function) or P (biological
process).

KEGG pathways are NOT taken from Swiss-Prot's KEGG cross-ref, which points at organism gene
IDs rather than pathways. Instead EC numbers are mapped to KEGG pathways through the KEGG
REST API (enzyme -> pathway), which is appropriate for a non-model organism with no KEGG
genome entry of its own.
"""
import argparse, collections, csv, gzip, pathlib, sys, time, urllib.request

ROOT = pathlib.Path(__file__).resolve().parents[1]
ASPECT = {"C": "GO_CC", "F": "GO_MF", "P": "GO_BP"}


def scan_swissprot(dat, wanted):
    info = collections.defaultdict(lambda: {"go": [], "ec": set()})
    acc, keep = None, False
    with gzip.open(dat, "rt", errors="replace") as fh:
        for line in fh:
            if line.startswith("AC "):
                if acc is None:
                    for a in line[5:].replace(" ", "").rstrip(";\n").split(";"):
                        if a in wanted:
                            acc, keep = a, True
                            break
                    else:
                        acc, keep = "?", False
            elif keep and line.startswith("DR   GO;"):
                parts = [p.strip() for p in line.split(";")]
                if len(parts) >= 3 and len(parts[2]) > 2 and parts[2][1] == ":":
                    asp, name = parts[2][0], parts[2][2:]
                    if asp in ASPECT:
                        info[acc]["go"].append((ASPECT[asp], parts[1], name))
            elif keep and line.startswith("DE ") and "EC=" in line:
                for frag in line.split("EC=")[1:]:
                    e = frag.split(";")[0].strip().rstrip(",")
                    if e and "-" not in e:
                        info[acc]["ec"].add(e)
            elif line.startswith("//"):
                acc, keep = None, False
    return info


def kegg_ec_to_pathway(ecs):
    """EC -> KEGG pathway via the KEGG REST link endpoint, plus pathway names."""
    ec2pw, pw_name = collections.defaultdict(set), {}
    try:
        with urllib.request.urlopen("https://rest.kegg.jp/list/pathway", timeout=90) as r:
            for line in r.read().decode().splitlines():
                p = line.split("\t")
                if len(p) == 2:
                    pw_name[p[0].replace("path:", "")] = p[1]
    except Exception as e:
        print(f"  WARNING: could not fetch KEGG pathway names ({e})")
    try:
        with urllib.request.urlopen("https://rest.kegg.jp/link/pathway/enzyme", timeout=180) as r:
            for line in r.read().decode().splitlines():
                p = line.split("\t")
                if len(p) == 2:
                    ec = p[0].replace("ec:", "")
                    # both KEGG endpoints return "map" ids; earlier code rewrote them to "ko"
                    # on one side only, which broke every name lookup
                    pw = p[1].replace("path:", "")
                    if ec in ecs:
                        ec2pw[ec].add(pw)
    except Exception as e:
        print(f"  WARNING: could not fetch EC->pathway links ({e})")
    return ec2pw, pw_name


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--functional", default=str(ROOT / "results/annotation/bom_ss5_functional.tsv"))
    ap.add_argument("--gtf", default=str(ROOT / "refs/BOM_ss5/BOM_ss5_genomic.gtf"))
    ap.add_argument("--dat", default=str(ROOT / "refs/uniprot/uniprot_sprot.dat.gz"))
    ap.add_argument("--out", default=str(ROOT / "results/annotation/bom_ss5_terms.tsv"))
    args = ap.parse_args()

    rows = list(csv.DictReader(open(args.functional), delimiter="\t"))
    prot2acc = {r["protein_id"]: r["sprot_acc"] for r in rows if r.get("sprot_acc")}
    print(f"proteins with Swiss-Prot hit: {len(prot2acc):,}")

    info = scan_swissprot(args.dat, set(prot2acc.values()))
    print(f"accessions resolved: {len(info):,}")

    all_ec = {e for v in info.values() for e in v["ec"]}
    print(f"distinct EC numbers: {len(all_ec):,}; querying KEGG ...")
    ec2pw, pw_name = kegg_ec_to_pathway(all_ec)
    print(f"EC numbers with a KEGG pathway: {len(ec2pw):,}; pathways named: {len(pw_name):,}")

    # gene_id <- protein_id from the stock GTF CDS records
    import re
    gene_of = {}
    for line in open(args.gtf):
        if line.startswith("#"): continue
        p = line.split("\t")
        if len(p) < 9 or p[2] != "CDS": continue
        g = re.search(r'gene_id "([^"]+)"', p[8]); pr = re.search(r'protein_id "([^"]+)"', p[8])
        if g and pr: gene_of.setdefault(pr.group(1), g.group(1))

    out = pathlib.Path(args.out)
    n = 0
    with out.open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["gene_id", "protein_id", "ontology", "term_id", "term_name"])
        for prot, acc in sorted(prot2acc.items()):
            gid = gene_of.get(prot)
            if not gid: continue
            rec = info.get(acc)
            if not rec: continue
            seen = set()
            for onto, tid, name in rec["go"]:
                if (onto, tid) in seen: continue
                seen.add((onto, tid)); w.writerow([gid, prot, onto, tid, name]); n += 1
            for e in rec["ec"]:
                for pw in ec2pw.get(e, ()):
                    if ("KEGG", pw) in seen: continue
                    seen.add(("KEGG", pw))
                    w.writerow([gid, prot, "KEGG", pw, pw_name.get(pw, pw)]); n += 1

    print(f"term assignments written: {n:,}")
    per = collections.Counter()
    for line in open(out):
        p = line.split("\t")
        if p[0] != "gene_id": per[p[2]] += 1
    print("  by ontology:", dict(per))
    print(f"wrote {out}")


if __name__ == "__main__":
    sys.exit(main())
