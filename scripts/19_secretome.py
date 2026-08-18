#!/usr/bin/env python3
"""Phase 5 -- secretome inference by homology transfer from Swiss-Prot.

De novo signal-peptide prediction was not available on this machine: SignalP and Phobius are
licence-restricted, and DeepSig's conda build pins tensorflow 2.2.0, which has no osx-arm64
package. Rather than roll a bespoke predictor, this transfers the CURATED SIGNAL and
TRANSMEM feature annotations from the best-hit Swiss-Prot entry of each protein.

A protein is called secreted when its Swiss-Prot homologue carries a SIGNAL feature and has
no TRANSMEM helix outside that signal region (the standard signal-yes/TM-no rule).

Limitations, stated plainly: this only covers proteins WITH a Swiss-Prot hit (~41% of the
proteome), it inherits the homologue's annotation rather than inspecting the fungal
sequence, and it will miss lineage-specific secreted proteins that have no characterised
homologue -- a real concern in basidiomycetes. Treat counts as a lower bound.
"""
import csv, gzip, pathlib, re, sys, collections

ROOT = pathlib.Path(__file__).resolve().parents[1]
func = sys.argv[1] if len(sys.argv) > 1 else str(ROOT / "results/annotation/bom_ss5_functional.tsv")
dat  = ROOT / "refs/uniprot/uniprot_sprot.dat.gz"
out  = pathlib.Path(sys.argv[2] if len(sys.argv) > 2 else ROOT / "results/annotation/bom_ss5_secretome.tsv")

rows = list(csv.DictReader(open(func), delimiter="\t"))
want = {r["sprot_acc"] for r in rows if r.get("sprot_acc")}
print(f"proteins with a Swiss-Prot hit: {len(rows):,}; accessions to scan: {len(want):,}")

feats = collections.defaultdict(lambda: {"signal": None, "tm": []})
acc, keep = None, False
with gzip.open(dat, "rt", errors="replace") as fh:
    for line in fh:
        if line.startswith("AC "):
            if acc is None:
                for a in line[5:].replace(" ", "").rstrip(";\n").split(";"):
                    if a in want:
                        acc, keep = a, True
                        break
                else:
                    acc, keep = "?", False
        elif keep and line.startswith("FT   SIGNAL"):
            m = re.search(r"(\d+)\.\.(\d+)", line)
            if m: feats[acc]["signal"] = (int(m.group(1)), int(m.group(2)))
        elif keep and line.startswith("FT   TRANSMEM"):
            m = re.search(r"(\d+)\.\.(\d+)", line)
            if m: feats[acc]["tm"].append((int(m.group(1)), int(m.group(2))))
        elif line.startswith("//"):
            acc, keep = None, False

n_sec = 0
with out.open("w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t")
    w.writerow(["protein_id", "sprot_acc", "pident", "protein_name", "EC",
                "signal_peptide", "n_TM_outside_signal", "secreted"])
    for r in rows:
        f = feats.get(r["sprot_acc"], {"signal": None, "tm": []})
        sig = f["signal"]
        tm_out = 0 if not sig else sum(1 for s, e in f["tm"] if s > sig[1])
        if sig is None:
            tm_out = len(f["tm"])
        secreted = bool(sig) and tm_out == 0
        n_sec += secreted
        w.writerow([r["protein_id"], r["sprot_acc"], r["pident"], r["protein_name"], r["EC"],
                    f"{sig[0]}-{sig[1]}" if sig else "", tm_out, int(secreted)])

print(f"predicted secreted (signal peptide, no TM beyond it): {n_sec:,} "
      f"({100*n_sec/len(rows):.1f}% of annotated proteins)")
print(f"wrote {out}")
