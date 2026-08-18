#!/usr/bin/env python3
"""Phase 0 -- build a GeneLab-style runsheet from Myco_Seq_Data_Plan.csv.

Every FASTQ path is verified to exist on disk. The sample sheet was corrected once already
(the well IDs for the two middle tissue groups were transposed), so this asserts rather than
trusts: an unmatched or duplicated filename is a hard error, not a warning.
"""
import csv
import pathlib
import sys

ROOT = pathlib.Path(__file__).resolve().parents[1]
PLAN = ROOT.parent / "Myco_Seq_Data_Plan.csv"
FASTQ_DIR = ROOT.parent / "GPNJ7M_fastq"
OUT = ROOT / "metadata" / "runsheet.csv"

rows = []
with PLAN.open() as fh:
    for rec in csv.DictReader(fh):
        tissue = (rec.get("Tissue") or "").strip()
        fname = (rec.get("File name") or "").strip()
        if not tissue or not fname:
            continue
        path = FASTQ_DIR / fname
        if not path.exists():
            sys.exit(f"ERROR: {fname} listed in the plan but absent from {FASTQ_DIR}")
        sample = fname.replace(".fastq.gz", "")
        rows.append(
            {
                "sample_name": sample,
                "well": sample.split("_sample_")[1],
                "index": int(sample.split("_")[1]),
                "Factor Value[Tissue]": tissue,
                "organism": "Pleurotus ostreatus",
                "paired_end": "false",
                "has_ERCC": "false",
                "read1_path": str(path),
            }
        )

on_disk = {p.name for p in FASTQ_DIR.glob("*.fastq.gz")}
listed = {r["read1_path"].rsplit("/", 1)[1] for r in rows}
if on_disk != listed:
    sys.exit(f"ERROR: plan/disk mismatch. only-on-disk={on_disk - listed} only-in-plan={listed - on_disk}")
if len(listed) != len(rows):
    sys.exit("ERROR: duplicate filenames in the plan")

rows.sort(key=lambda r: r["index"])
OUT.parent.mkdir(parents=True, exist_ok=True)
with OUT.open("w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0]))
    w.writeheader()
    w.writerows(rows)

print(f"wrote {OUT} ({len(rows)} samples, all paths verified)")
counts = {}
for r in rows:
    counts.setdefault(r["Factor Value[Tissue]"], []).append(r["well"])
for tissue, wells in counts.items():
    print(f"  {tissue:<20} n={len(wells)}  wells={','.join(wells)}")
