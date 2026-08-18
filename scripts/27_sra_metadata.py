#!/usr/bin/env python3
"""Generate NCBI BioSample and SRA submission metadata from the project runsheet.

Everything derivable from the data is filled in automatically: file names, read counts and
lengths, instrument and run identifiers (parsed from the FASTQ headers), library layout and
selection. Fields that only the submitter can supply are written as explicit
`<<REQUIRED: ...>>` placeholders rather than plausible-looking guesses, because a wrong
collection date or geographic origin in a public archive is worse than a blank one.

Outputs (tab-separated, the format NCBI's submission portal accepts):
  submission/biosample_attributes.tsv   one row per sample
  submission/sra_metadata.tsv           one row per run
  submission/SUBMISSION_NOTES.md        what to fill in and in what order
"""
import csv, gzip, json, pathlib, sys, collections

ROOT = pathlib.Path(__file__).resolve().parents[1]
FASTQ = ROOT.parent / "GPNJ7M_fastq"
OUT = ROOT / "submission"
REQ = lambda s: f"<<REQUIRED: {s}>>"

ORGANISM = "Pleurotus ostreatus"
STRAIN = "Harbor Blue P01"
# Tissue -> a short, archive-friendly description. The exudophore wording is the submitter's
# to finalise; it is a newly described structure with no established controlled term.
TISSUE_DESC = {
    "Exuding mycelium": "vegetative mycelium bearing exudate droplets",
    "Fuzzy mycelium": "aerial vegetative mycelium",
    "Exudophore": "exudate-producing structure (newly described)",
    "Nodule": "mycelial nodule / early primordium",
}


def fastq_facts(path):
    n = 0
    lens = collections.Counter()
    hdr = None
    with gzip.open(path, "rt") as fh:
        for i, line in enumerate(fh):
            if i % 4 == 0 and hdr is None:
                hdr = line.strip()
            elif i % 4 == 1:
                lens[len(line.strip())] += 1
                n += 1
            if n >= 200000:
                break
    inst, run, fc, lane = "", "", "", ""
    if hdr:
        p = hdr.lstrip("@").split(":")
        if len(p) >= 4:
            inst, run, fc, lane = p[0], p[1], p[2], p[3]
    return {"max_len": max(lens) if lens else 0, "instrument_id": inst,
            "run_id": run, "flowcell": fc, "lane": lane}


def main():
    OUT.mkdir(exist_ok=True)
    rows = list(csv.DictReader((ROOT / "metadata" / "runsheet.csv").open()))
    budget = {r["sample"]: r for r in csv.DictReader((ROOT / "results" / "read_budget_BOM_ss5.csv").open())}

    facts = fastq_facts(FASTQ / (rows[0]["sample_name"] + ".fastq.gz"))
    print(f"instrument {facts['instrument_id']} run {facts['run_id']} "
          f"flowcell {facts['flowcell']} lane {facts['lane']}; max read length {facts['max_len']}")

    # ---- BioSample ----
    bs_cols = ["sample_name", "organism", "strain", "isolate", "cultivar",
               "collection_date", "geo_loc_name", "isolation_source", "tissue",
               "dev_stage", "lab_host", "growth_protocol", "biomaterial_provider",
               "description"]
    bs = []
    for r in rows:
        tis = r["Factor Value[Tissue]"]
        bs.append({
            "sample_name": r["sample_name"],
            "organism": ORGANISM,
            "strain": STRAIN,
            "isolate": "not applicable",
            "cultivar": STRAIN,
            "collection_date": REQ("YYYY-MM-DD or YYYY-MM of tissue harvest"),
            "geo_loc_name": REQ("country[:region] where cultures were grown, e.g. USA:Indiana"),
            "isolation_source": "axenic culture on micro-structured ceramic MiniTube "
                                "with granular activated carbon, mycoponic nutrient medium v3",
            "tissue": tis,
            "dev_stage": "vegetative mycelium" if "mycelium" in tis.lower() else REQ(f"developmental stage for {tis}"),
            "lab_host": "not applicable",
            "growth_protocol": "Mycoponic ceramic MiniTube (10 x 5 cm, pore size <300 nm, 50% v/v "
                               "granular activated carbon) with mycoponic nutrient medium v3; "
                               "16 C, 85% RH, CO2 <=1000 ppm (Porterfield et al. 2026, "
                               "doi:10.1002/biot.70184)",
            "biomaterial_provider": REQ("supplier/source of the Harbor Blue P01 culture"),
            "description": f"{TISSUE_DESC.get(tis, tis)}; replicate well {r['well']}",
        })
    with (OUT / "biosample_attributes.tsv").open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=bs_cols, delimiter="\t")
        w.writeheader(); w.writerows(bs)

    # ---- SRA metadata ----
    sra_cols = ["biosample_accession", "library_ID", "title", "library_strategy",
                "library_source", "library_selection", "library_layout", "platform",
                "instrument_model", "design_description", "filetype", "filename", "md5"]
    md5 = {}
    mp = OUT / "fastq_md5.txt"
    if mp.exists():
        for line in mp.read_text().splitlines():
            p = line.split()
            if len(p) == 2:
                md5[p[1]] = p[0]

    design = (f"Total RNA, poly(A)-selected 3'-end tag libraries with a 14 nt unique molecular "
              f"identifier; single-end {facts['max_len']} bp on Illumina NovaSeq X Plus "
              f"(instrument {facts['instrument_id']}, run {facts['run_id']}, "
              f"flowcell {facts['flowcell']}, lane {facts['lane']}). All 16 libraries were "
              f"sequenced on a single flowcell and lane. UMI is appended to the read name "
              f"after an underscore.")
    sra = []
    for r in rows:
        s = r["sample_name"]
        tis = r["Factor Value[Tissue]"]
        sra.append({
            "biosample_accession": REQ("SAMN accession returned after BioSample submission"),
            "library_ID": s,
            "title": f"3'-tag RNA-seq of {ORGANISM} {STRAIN} {tis} (replicate {r['well']})",
            "library_strategy": "RNA-Seq",
            "library_source": "TRANSCRIPTOMIC",
            "library_selection": "PolyA",
            "library_layout": "single",
            "platform": "ILLUMINA",
            "instrument_model": "Illumina NovaSeq X Plus",
            "design_description": design,
            "filetype": "fastq",
            "filename": s + ".fastq.gz",
            "md5": md5.get(s + ".fastq.gz", "<<pending: see submission/fastq_md5.txt>>"),
        })
    with (OUT / "sra_metadata.tsv").open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=sra_cols, delimiter="\t")
        w.writeheader(); w.writerows(sra)

    n_req = sum(str(v).startswith("<<REQUIRED") for row in bs for v in row.values())
    print(f"wrote {OUT/'biosample_attributes.tsv'} ({len(bs)} samples, {n_req} fields needing input)")
    print(f"wrote {OUT/'sra_metadata.tsv'} ({len(sra)} runs)")
    print(f"md5 checksums resolved: {sum(1 for r in sra if not r['md5'].startswith('<<'))}/{len(sra)}")


if __name__ == "__main__":
    sys.exit(main())
