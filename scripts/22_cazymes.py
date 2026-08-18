#!/usr/bin/env python3
"""Phase 5 -- CAZyme-family assignment from Pfam domains.

dbCAN is the standard tool, but its download service (bcb.unl.edu/dbCAN2/download) now
redirects to a landing page and no database file was retrievable. Pfam-A covers the CAZy
families directly, so domains are assigned with HMMER against Pfam-A and mapped to CAZy
classes by curated keyword rules over the Pfam descriptions.

This is labelled honestly throughout as **Pfam-derived**, not dbCAN. It will differ from a
dbCAN run: dbCAN uses CAZy-specific HMMs with family-specific thresholds, whereas this uses
Pfam families and a global E-value. Family-level calls (e.g. "GH5") are therefore not
attempted -- only CAZy CLASS (GH/GT/CE/PL/AA/CBM), which the Pfam descriptions support.
"""
import argparse, collections, csv, gzip, pathlib, re, sys

ROOT = pathlib.Path(__file__).resolve().parents[1]

# CAZy class <- keyword rules over Pfam DESC/NAME. Ordered: first match wins.
RULES = [
    ("AA",  r"lytic polysaccharide monooxygenase|LPMO|laccase|multicopper oxidase|"
            r"lignin peroxidase|manganese peroxidase|versatile peroxidase|"
            r"dye-decolou?ri[sz]ing peroxidase|DyP|glyoxal oxidase|GMC oxidoreductase|"
            r"cellobiose dehydrogenase|aryl-alcohol oxidase|copper radical oxidase|"
            r"galactose oxidase|vanillyl-alcohol oxidase|benzoquinone reductase"),
    ("GH",  r"glycosyl(?: |-)?hydrolase|glycoside hydrolase|cellulase|chitinase|"
            r"xylanase|glucanase|amylase|mannanase|galactosidase|glucosidase|"
            r"chitosanase|lysozyme|invertase|trehalase|pectinase|arabinofuranosidase"),
    ("GT",  r"glycosyl(?: |-)?transferase|glycosyltransferase"),
    ("PL",  r"polysaccharide lyase|pectate lyase|alginate lyase|rhamnogalacturonan lyase|"
            r"chondroitin lyase|heparinase"),
    ("CE",  r"carbohydrate esterase|acetyl xylan esterase|pectin methylesterase|"
            r"chitin deacetylase|rhamnogalacturonan acetylesterase|cutinase"),
    ("CBM", r"carbohydrate(?: |-)binding module|cellulose(?: |-)binding|chitin(?: |-)binding|"
            r"starch(?: |-)binding|CBM_|carbohydrate binding"),
]
COMPILED = [(c, re.compile(p, re.I)) for c, p in RULES]


def pfam_descriptions(hmm_path):
    """NAME/ACC/DESC triples straight from the HMM library."""
    desc = {}
    name = acc = None
    op = gzip.open if str(hmm_path).endswith(".gz") else open
    with op(hmm_path, "rt", errors="replace") as fh:
        for line in fh:
            if line.startswith("NAME "): name = line[5:].strip()
            elif line.startswith("ACC  "): acc = line[5:].strip()
            elif line.startswith("DESC "):
                d = line[5:].strip()
                if name: desc[name] = (acc or "", d)
            elif line.startswith("//"): name = acc = None
    return desc


def classify(name, d):
    text = f"{name} {d}"
    for cls, rx in COMPILED:
        if rx.search(text):
            return cls
    return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tbl", default=str(ROOT / "results/annotation/bom_ss5_pfam.tbl"))
    ap.add_argument("--hmm", default=str(ROOT / "refs/pfam/Pfam-A.hmm"))
    ap.add_argument("--out", default=str(ROOT / "results/annotation/bom_ss5_cazymes.tsv"))
    args = ap.parse_args()

    desc = pfam_descriptions(args.hmm)
    print(f"Pfam models with descriptions: {len(desc):,}")
    cazy_fams = {n: (a, d, classify(n, d)) for n, (a, d) in desc.items() if classify(n, d)}
    print(f"Pfam families matching a CAZy class: {len(cazy_fams):,}")
    print("  by class:", dict(collections.Counter(v[2] for v in cazy_fams.values())))

    best = {}
    with open(args.tbl) as fh:
        for line in fh:
            if line.startswith("#"): continue
            p = line.split()
            if len(p) < 6: continue
            prot, fam, ev = p[0], p[2], float(p[4])
            if fam not in cazy_fams: continue
            if prot not in best or ev < best[prot][1]:
                best[prot] = (fam, ev)

    out = pathlib.Path(args.out)
    with out.open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["protein_id", "pfam_name", "pfam_acc", "cazy_class", "evalue", "description"])
        for prot, (fam, ev) in sorted(best.items()):
            acc, d, cls = cazy_fams[fam]
            w.writerow([prot, fam, acc, cls, f"{ev:.1e}", d])

    print(f"\nproteins with a CAZy-class domain: {len(best):,}")
    print("  by class:", dict(collections.Counter(cazy_fams[f][2] for f, _ in best.values())))
    print(f"wrote {out}")


if __name__ == "__main__":
    sys.exit(main())
