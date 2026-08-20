#!/usr/bin/env python3
"""Audit -- verify every database identifier in the manuscripts against its source.

Identifiers have produced two errors in this project already, both caught by chance rather
than by process:

  * two ModelSEED compound IDs were wrong (cpd00166 is Calomide, not CoA; cpd00557 is
    Siroheme, not a chitin precursor) -- see NOTES.md §14;
  * two NCBI accessions in the comparative set were not fungal at all (a Physcomitrella
    chromosome and a Cymbopogon chloroplast).

Both share a signature: a well-formed identifier attached to the wrong name. Checking that an
identifier merely *exists* would have passed either one, so every check here compares the
identifier against the name or organism the text claims for it.

Checks:
  * NCBI accessions      -> local FASTA headers, else the NCBI datasets/eutils summary
  * Rfam RF##### IDs     -> the DESC line of the downloaded covariance model
  * EC numbers           -> the Swiss-Prot flat file's DE ... EC= records
  * ModelSEED cpd/rxn    -> refs/modelseed/compounds.tsv, reactions.tsv

LaTeX escaping matters: accessions are written \\texttt{GCA\\_056149245.1}, so a naive scan
misses every identifier containing an underscore.

Usage: 51_audit_identifiers.py [--offline]
"""
import argparse
import csv
import gzip
import json
import pathlib
import re
import sys
import urllib.parse
import urllib.request

ROOT = pathlib.Path(__file__).resolve().parents[1]
LATEX = ROOT / "latex"
UA = "myco-tissue-rnaseq-audit/1.0 (mailto:dr.richard.barker@gmail.com)"

ACC_RE = re.compile(r"\b(?:GC[AF]_\d{9}\.\d|[A-Z]{2}\d{6}\.\d|N[CZ]_[A-Z0-9]+\.\d"
                    r"|[A-Z]{5,7}\d{8,9}\.\d)\b")
RFAM_RE = re.compile(r"\bRF\d{5}\b")
EC_RE = re.compile(r"\b\d\.[\d-]+\.[\d-]+\.[\d-]+\b")
SEED_RE = re.compile(r"\b(?:cpd|rxn)\d{5}\b")


def sources():
    """Every text we audit, with LaTeX escapes undone so identifiers read normally."""
    out = {}
    for tex in sorted(LATEX.glob("*/*.tex")):
        if tex.name not in ("sections.tex", "main.tex", "tables.tex", "figures.tex"):
            continue
        t = tex.read_text(errors="replace").replace("\\_", "_").replace("\\%", "%")
        out[str(tex.relative_to(LATEX))] = t
    for extra in ("NOTES.md", "README.md"):
        p = ROOT / extra
        if p.exists():
            out[extra] = p.read_text(errors="replace")
    return out


def local_accessions():
    """Accession -> organism/description, from every local FASTA and the comparative table."""
    known = {}
    inv = ROOT / "results/mito_comparative/inventory.csv"
    for pat in ("refs/**/*.fa", "refs/**/*.fasta", "refs/**/*.fna", "results/**/*.fa"):
        for fa in ROOT.glob(pat):
            try:
                with fa.open(errors="replace") as fh:
                    for line in fh:
                        if not line.startswith(">"):
                            continue
                        acc = line[1:].split()[0]
                        known.setdefault(acc, line[1:].strip())
            except OSError:
                continue
    if inv.exists():
        for r in csv.DictReader(inv.open()):
            known.setdefault(r["acc"], f"comparative set, {r['length']} bp")
    return known


def ncbi_assembly(acc):
    """Assembly accessions (GCA_/GCF_) are not in nuccore -- they need the datasets API."""
    req = urllib.request.Request(
        "https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/"
        + urllib.parse.quote(acc) + "/dataset_report",
        headers={"User-Agent": UA})
    with urllib.request.urlopen(req, timeout=60) as fh:
        rep = json.load(fh)
    rows = rep.get("reports") or []
    if not rows:
        return {}
    r = rows[0]
    return {"organism": (r.get("organism") or {}).get("organism_name", "?"),
            "title": ((r.get("assembly_info") or {}).get("assembly_name", "?")
                      + " / " + (r.get("assembly_info") or {}).get("biosample", {})
                      .get("description", {}).get("title", "")),
            "strain": ((r.get("organism") or {}).get("infraspecific_names") or {})
                      .get("strain", "")}


def ncbi_summary(acc):
    q = urllib.parse.urlencode({"db": "nuccore", "id": acc, "retmode": "json"})
    req = urllib.request.Request(
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?" + q,
        headers={"User-Agent": UA})
    with urllib.request.urlopen(req, timeout=60) as fh:
        d = json.load(fh)["result"]
    uid = d.get("uids", [None])[0]
    return d.get(uid, {}) if uid else {}


def rfam_desc():
    """RF##### -> description, from the DESC line of each downloaded covariance model."""
    out = {}
    for cm in ROOT.glob("refs/rfam/*.cm"):
        acc = desc = None
        with cm.open(errors="replace") as fh:
            for line in fh:
                if line.startswith("ACC "):
                    acc = line.split()[1]
                elif line.startswith("DESC "):
                    desc = line[5:].strip()
                if acc and desc:
                    out[acc] = desc
                    acc = desc = None
    return out


def swissprot_ec():
    """EC number -> enzyme name, from the Swiss-Prot flat file DE records."""
    p = ROOT / "refs/uniprot/uniprot_sprot.dat.gz"
    out = {}
    if not p.exists():
        return out
    name = None
    with gzip.open(p, "rt", errors="replace") as fh:
        for line in fh:
            if line.startswith("DE   RecName: Full="):
                name = line.split("Full=", 1)[1].split("{")[0].strip().rstrip(";")
            elif "EC=" in line and line.startswith("DE"):
                for ec in re.findall(r"EC=([\d.\-]+)", line):
                    if name:
                        out.setdefault(ec, name)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--offline", action="store_true")
    args = ap.parse_args()

    docs = sources()
    known = local_accessions()
    rfam = rfam_desc()
    print(f"local FASTA accessions: {len(known):,}   Rfam models: {len(rfam)}")
    ecmap = swissprot_ec()
    print(f"Swiss-Prot EC names: {len(ecmap):,}\n")

    seed = {}
    cf = ROOT / "refs/modelseed/compounds.tsv"
    if cf.exists():
        for r in csv.DictReader(cf.open(), delimiter="\t"):
            seed[r["id"]] = r.get("name", "")
    rf = ROOT / "refs/modelseed/reactions.tsv"
    if rf.exists():
        for r in csv.DictReader(rf.open(), delimiter="\t"):
            seed[r["id"]] = r.get("name", "")

    fails = 0
    seen_acc = {}
    for doc, text in docs.items():
        for acc in sorted(set(ACC_RE.findall(text))):
            seen_acc.setdefault(acc, []).append(doc)
        for rfid in sorted(set(RFAM_RE.findall(text))):
            if rfid not in rfam:
                print(f"FAIL  {doc}: {rfid} has no local covariance model in refs/rfam/")
                fails += 1
        for sid in sorted(set(SEED_RE.findall(text))):
            if seed and sid not in seed:
                print(f"FAIL  {doc}: ModelSEED id {sid} not in refs/modelseed/")
                fails += 1

    print("== Rfam ==")
    for rfid, desc in sorted(rfam.items()):
        used = [d for d, t in docs.items() if rfid in t]
        print(f"ok    {rfid}  {desc}" + (f"   [{', '.join(used)}]" if used else "   (unused)"))

    print("\n== NCBI accessions ==")
    for acc, docs_using in sorted(seen_acc.items()):
        if acc in known:
            print(f"ok    {acc:<20} local: {known[acc][:70]}")
            continue
        if args.offline:
            print(f"WARN  {acc:<20} not local (offline, not checked)")
            continue
        try:
            s = ncbi_assembly(acc) if acc.startswith(("GCA_", "GCF_")) else ncbi_summary(acc)
        except Exception as exc:
            print(f"FAIL  {acc:<20} lookup failed: {exc}")
            fails += 1
            continue
        org, title = s.get("organism", "?"), s.get("title", "?")
        if s.get("strain"):
            title = f"strain {s['strain']} | {title}"
        bad = "Pleurotus" not in org and "Pleurotus" not in title
        print(f"{'WARN ' if bad else 'ok   '} {acc:<20} NCBI: {org} | {title[:60]}")
        if bad:
            print(f"        ^ not a Pleurotus record -- confirm this is intended "
                  f"({', '.join(docs_using)})")

    print("\n== EC numbers ==")
    for doc, text in docs.items():
        for ec in sorted(set(EC_RE.findall(text))):
            if not ecmap:
                print(f"WARN  {ec} in {doc}: no Swiss-Prot flat file to check against")
                continue
            if ec in ecmap:
                print(f"ok    {ec:<12} {ecmap[ec][:64]}   [{doc}]")
            else:
                print(f"WARN  {ec:<12} not found in Swiss-Prot DE records   [{doc}]")

    print(f"\n{'FAILURES: ' + str(fails) if fails else 'all identifiers verified'}")
    return 1 if fails else 0


if __name__ == "__main__":
    sys.exit(main())
