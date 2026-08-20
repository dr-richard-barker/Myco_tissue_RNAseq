#!/usr/bin/env python3
"""Audit -- verify every bibliography entry against CrossRef, and every \\cite against the bib.

Written after an entry in refs.bib was found to be fabricated: DOI 10.3389/fmicb.2020.00863
was recorded as "Zaccaron & Bluhm, Erosion of Genome Sequences and Emergence of Novel Genes..."
when it is in fact Song, Geng & Li, "The Mitochondrial Genome of the Phytopathogenic Fungus
Bipolaris sorokiniana...". The DOI was right and everything attached to it was wrong, so a
check that merely resolved the DOI would have passed it. This compares the metadata.

Checks, per entry:
  * DOI resolves at CrossRef
  * title matches (token overlap; punctuation, case and LaTeX accents ignored)
  * every author surname in the entry appears in CrossRef's author list, and vice versa
  * year matches within 1 (online-first vs issue year), volume and first page match

Plus, across the LaTeX sources:
  * every \\cite{key} resolves to an entry in refs.bib
  * every entry in refs.bib is cited at least once

Entries with no DOI are listed, not failed: barrnap and the GeneLab pipeline spec genuinely
have none.

Usage:  50_audit_citations.py [--bib latex/refs.bib] [--offline]
Exit status is non-zero if any check fails, so it can gate the build.
"""
import argparse
import json
import pathlib
import re
import sys
import unicodedata
import time
import urllib.parse
import urllib.request

ROOT = pathlib.Path(__file__).resolve().parents[1]
UA = "myco-tissue-rnaseq-audit/1.0 (mailto:dr.richard.barker@gmail.com)"
NO_DOI_OK = {"seemann2013barrnap", "genelab7101G"}

ACCENTS = {r"\'e": "e", r"\'a": "a", r"\'o": "o", r"\'i": "i", r"\'u": "u",
           r'\"o': "o", r'\"u': "u", r'\"a': "a", r"\~n": "n", r"\c c": "c"}


def deaccent(s):
    """LaTeX accent macros and Unicode combining marks both flattened to ASCII, so that
    F{\\'e}randon in the bib compares equal to Ferandon in CrossRef's Unicode."""
    for k, v in ACCENTS.items():
        s = s.replace("{" + k + "}", v).replace(k, v)
    s = unicodedata.normalize("NFKD", s)
    return "".join(c for c in s if not unicodedata.combining(c))


def norm_title(s):
    s = deaccent(s).lower()
    s = re.sub(r"<[^>]+>", " ", s)
    s = re.sub(r"[^a-z0-9 ]", " ", s)
    return [w for w in s.split() if w]


def parse_bib(path):
    """Minimal bibtex reader. Brace-matched rather than line-based: entries in this file end
    with the last field and the entry brace on one line ("...year={2020}}"), which a regex
    anchored on a newline silently skips."""
    text = path.read_text()
    entries = []
    for m in re.finditer(r"@(\w+)\s*\{\s*([^,\s]+)\s*,", text):
        depth, i = 1, m.end()
        while i < len(text) and depth:
            if text[i] == "{":
                depth += 1
            elif text[i] == "}":
                depth -= 1
            i += 1
        body = text[m.end():i - 1]
        fields = {}
        for fm in re.finditer(r"(\w+)\s*=\s*\{((?:[^{}]|\{[^{}]*\})*)\}", body):
            fields[fm.group(1).lower()] = fm.group(2).strip()
        entries.append({"key": m.group(2).strip(), "type": m.group(1), **fields})
    return entries


def crossref(doi):
    url = "https://api.crossref.org/works/" + urllib.parse.quote(doi)
    req = urllib.request.Request(url, headers={"User-Agent": UA})
    with urllib.request.urlopen(req, timeout=60) as fh:
        return json.load(fh)["message"]


def surnames(bib_author):
    """Surnames from a bibtex author field, plus flags for the two idioms that are not names.

    Returns (surnames, truncated, corporate) where `truncated` means the field ended in
    "and others" (bibtex for et al., so CrossRef listing more authors is expected) and
    `corporate` collects braced institutional authors such as {UniProt Consortium}.
    """
    out, truncated, corporate = [], False, []
    for part in re.split(r"\s+and\s+", bib_author or ""):
        part = deaccent(part).strip()
        if not part:
            continue
        if part.lower() == "others":
            truncated = True
            continue
        if part.startswith("{") and part.endswith("}"):
            corporate.append(part.strip("{}").lower())
            continue
        out.append((part.split(",")[0] if "," in part else part.split()[-1]).lower())
    return out, truncated, corporate


def check_entry(e, meta):
    """Return (problems, warnings) for one entry."""
    bad, warn = [], []

    got = norm_title(meta.get("title", [""])[0] if meta.get("title") else "")
    want = norm_title(e.get("title", ""))
    if want and got:
        overlap = len(set(want) & set(got)) / max(len(set(want)), 1)
        if overlap < 0.7:
            bad.append(f"TITLE  bib={' '.join(want)[:70]!r}\n            doi={' '.join(got)[:70]!r}")

    cr = [deaccent((a.get("family") or "")).lower()
          for a in meta.get("author", []) if a.get("family")]
    bib, truncated, corporate = surnames(e.get("author", ""))
    if cr and bib:
        # The first author is the strongest signal and CrossRef always carries it: a
        # mismatch there is the signature of an entry attached to the wrong paper.
        if not (bib[0] in cr[0] or cr[0] in bib[0]):
            bad.append(f"AUTHOR first author bib={bib[0]!r} doi={cr[0]!r}  <- likely wrong paper")
        missing = [a for a in bib if not any(a in c or c in a for c in cr)]
        if missing:
            # CrossRef records for older papers are often truncated to the first author.
            # Only treat a bib author absent from CrossRef as an error when CrossRef
            # actually lists a full author list to be absent from.
            if len(cr) >= len(bib):
                bad.append(f"AUTHOR not in CrossRef record: {missing}")
            else:
                warn.append(f"AUTHOR CrossRef record lists only {len(cr)}; unverified: {missing}")
        if not truncated:
            extra = [c for c in cr if not any(c in b or b in c for b in bib)]
            if extra:
                bad.append(f"AUTHOR in CrossRef but not in bib: {extra[:8]}"
                           f"{' ...' if len(extra) > 8 else ''}")
    elif corporate:
        warn.append(f"AUTHOR corporate author {corporate}, not checked against CrossRef")

    parts = (meta.get("issued") or {}).get("date-parts") or [[None]]
    cr_year = parts[0][0]
    try:
        bib_year = int(e.get("year", "0"))
    except ValueError:
        bib_year = 0
    if cr_year and bib_year and abs(cr_year - bib_year) > 1:
        bad.append(f"YEAR   bib={bib_year} doi={cr_year}")

    if e.get("volume") and meta.get("volume") and e["volume"] != str(meta["volume"]):
        bad.append(f"VOLUME bib={e['volume']} doi={meta['volume']}")

    if e.get("pages") and meta.get("page"):
        fp = lambda s: re.split(r"[-–]+", s.replace("--", "-"))[0].strip()
        if fp(e["pages"]) != fp(str(meta["page"])):
            bad.append(f"PAGES  bib={e['pages']} doi={meta['page']}")
    return bad, warn


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bib", default=str(ROOT / "latex/refs.bib"))
    ap.add_argument("--texdir", default=str(ROOT / "latex"))
    ap.add_argument("--offline", action="store_true", help="skip CrossRef, check \\cite wiring only")
    args = ap.parse_args()

    entries = parse_bib(pathlib.Path(args.bib))
    print(f"refs.bib: {len(entries)} entries\n")
    failures = 0

    # ---- cite wiring -------------------------------------------------------
    keys = {e["key"] for e in entries}
    cited = set()
    OURS = ("sections.tex", "main.tex", "preamble.tex", "figures.tex", "tables.tex")
    for tex in sorted(pathlib.Path(args.texdir).rglob("*.tex")):
        if tex.name not in OURS:
            continue   # skip the vendored Springer class's own .tex files
        for m in re.finditer(r"\\cite[tp]?\{([^}]*)\}", tex.read_text(errors="replace")):
            cited.update(k.strip() for k in m.group(1).split(",") if k.strip())
    dangling = sorted(cited - keys)
    unused = sorted(keys - cited)
    if dangling:
        print(f"FAIL  {len(dangling)} \\cite key(s) with no bib entry: {dangling}")
        failures += len(dangling)
    else:
        print("ok    every \\cite key resolves to a bib entry")
    print(f"      {len(unused)} uncited entr(ies){': ' + str(unused) if unused else ''}\n")

    if args.offline:
        return 1 if failures else 0

    # ---- CrossRef ----------------------------------------------------------
    nodoi = []
    for e in entries:
        doi = e.get("doi")
        if not doi:
            nodoi.append(e["key"])
            continue
        try:
            meta = crossref(doi)
        except Exception as exc:
            print(f"FAIL  {e['key']:<26} DOI does not resolve ({doi}): {exc}")
            failures += 1
            time.sleep(0.3)
            continue
        bad, warn = check_entry(e, meta)
        if bad:
            failures += len(bad)
            print(f"FAIL  {e['key']:<26} {doi}")
        else:
            print(f"ok    {e['key']:<26} {doi}")
        for b in bad:
            print(f"        {b}")
        for w in warn:
            print(f"        note: {w}")
        time.sleep(0.3)

    for k in nodoi:
        flag = "ok   " if k in NO_DOI_OK else "WARN "
        print(f"{flag} {k:<26} no DOI" + ("" if k in NO_DOI_OK else "  <- unexpected"))
        if k not in NO_DOI_OK:
            failures += 1

    print(f"\n{'FAILURES: ' + str(failures) if failures else 'all citation checks passed'}")
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
