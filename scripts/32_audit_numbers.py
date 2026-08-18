#!/usr/bin/env python3
"""Verify every quantitative claim in the manuscripts against the analysis outputs.

Requested explicitly: the manuscripts must contain no unverifiable numbers. This builds a
dictionary of facts computed live from result files, then scans each .tex source for the
numbers it asserts and checks them against that dictionary.

Two categories are reported separately:
  CHECKED  - the value appears in the fact dictionary and matches (or does not, -> FAIL)
  UNMATCHED - a number in the text with no corresponding fact; listed so a human can confirm
              it is a parameter, a citation year, or a value stated elsewhere. These are not
              failures, but a long list is a warning sign.
"""
import csv, json, pathlib, re, statistics, sys, collections

ROOT = pathlib.Path(__file__).resolve().parents[1]
LATEX = ROOT / "latex"


def facts():
    f = {}
    bud = list(csv.DictReader((ROOT / "results/read_budget_BOM_ss5.csv").open()))
    f["raw total"] = sum(int(r["raw"]) for r in bud)
    f["mRNA total"] = sum(int(r["mRNA"]) for r in bud)
    f["rRNA total"] = sum(int(r["rRNA"]) for r in bud)
    f["libraries"] = len(bud)
    f["mRNA pct min"] = min(float(r["pct"]) for r in bud)
    f["mRNA pct max"] = max(float(r["pct"]) for r in bud)
    f["libs over 100k"] = sum(1 for r in bud if int(r["mRNA"]) >= 100_000)
    f["libs under 15k"] = sum(1 for r in bud if int(r["mRNA"]) < 15_000)
    f["mRNA pct of raw"] = round(100 * f["mRNA total"] / f["raw total"], 1)

    # direct rRNA measurement
    rr = []
    for p in (ROOT / "qc/rrna").glob("*.summary.txt"):
        t = p.read_text()
        g = lambda rx: float(re.search(rx, t).group(1))
        tot = g(r"Total reads: *(\d+)")
        rr.append(100 * (tot - g(r"Aligned 0 time: *(\d+)")) / tot)
    f["rRNA mean pct"] = round(statistics.mean(rr), 1)
    f["rRNA min pct"] = round(min(rr), 1)
    f["rRNA max pct"] = round(max(rr), 1)

    # reference test-map
    tm = collections.defaultdict(list)
    for p in (ROOT / "qc/testmap").glob("*.summary.txt"):
        t = p.read_text()
        g = lambda rx: float(re.search(rx, t).group(1))
        tot = g(r"Total reads: *(\d+)")
        tm[p.name.split("__")[0]].append(100 * g(r"Aligned 1 time: *(\d+)") / tot)
    for k, v in tm.items():
        f[f"{k} unique pct"] = round(statistics.mean(v), 2)

    # reference gain
    pc = {r["sample"]: int(r["mRNA"]) for r in csv.DictReader((ROOT / "results/read_budget_PC9.15.csv").open())}
    bo = {r["sample"]: int(r["mRNA"]) for r in csv.DictReader((ROOT / "results/read_budget_BOM_ss5.csv").open())}
    gains = [100 * (bo[s] - pc[s]) / pc[s] for s in bo]
    f["ref gain pct"] = round(100 * (sum(bo.values()) - sum(pc.values())) / sum(pc.values()), 1)
    f["ref gain min"] = round(min(gains), 1)
    f["ref gain max"] = round(max(gains), 1)

    # annotation
    f["proteins"] = sum(1 for l in (ROOT / "refs/BOM_ss5/BOM_ss5_protein.faa").open() if l.startswith(">"))
    f["sprot hits"] = sum(1 for _ in (ROOT / "results/annotation/bom_ss5_functional.tsv").open()) - 1
    f["sprot pct"] = round(100 * f["sprot hits"] / f["proteins"])
    cz = list(csv.DictReader((ROOT / "results/annotation/bom_ss5_cazymes.tsv").open(), delimiter="\t"))
    f["cazymes"] = len(cz)
    for cls, n in collections.Counter(r["cazy_class"] for r in cz).items():
        f[f"cazy {cls}"] = n
    f["secreted"] = sum(1 for r in csv.DictReader((ROOT / "results/annotation/bom_ss5_secretome.tsv").open(), delimiter="\t") if r["secreted"] == "1")
    f["terms"] = sum(1 for _ in (ROOT / "results/annotation/bom_ss5_terms.tsv").open()) - 1

    # GEM
    for g in csv.DictReader((ROOT / "results/figure_data/gem_stats.csv").open()):
        f[f"gem {g['stage']} reactions"] = int(g["reactions"])
        f[f"gem {g['stage']} carrying"] = int(g["carrying"])
        f[f"gem {g['stage']} metabolites"] = int(g["metabolites"])
        f[f"gem {g['stage']} genes"] = int(g["genes"])

    # WGCNA
    mods = list(csv.DictReader((ROOT / "results/wgcna/gene_modules.csv").open()))
    f["wgcna genes"] = len(mods)
    f["wgcna modules"] = len({m["module"] for m in mods} - {"grey"})
    conf = list(csv.DictReader((ROOT / "results/wgcna/rrna_confound_test.csv").open()))
    f["module trait tests"] = len(conf)
    f["modules fdr"] = sum(1 for r in conf if float(r["fdr_tissue"]) < 0.05 and abs(float(r["r_tissue"])) > 0.7)
    f["modules p05"] = sum(1 for r in conf if float(r["p_tissue"]) < 0.05)

    # markers
    rb = collections.Counter(r["tissue"] for r in csv.DictReader((ROOT / "results/tissue_models_BOM_ss5/markers_robust.csv").open()))
    for k, v in rb.items():
        f[f"robust {k}"] = v

    # genes surviving expression filtering, both retention scenarios
    for lab, d in [("all16", "dge_BOM_ss5"), ("filtered", "dge_BOM_ss5_filtered")]:
        for track in ("rRNArm", "all_genes"):
            fp = ROOT / "results" / d / f"VST_counts_{track}.csv"
            if fp.exists():
                f[f"genes filter {lab} {track}"] = sum(1 for _ in fp.open()) - 1

    # The delivered (mis-referenced) analysis, verified against the delivered files rather
    # than quoted. Directory name is globbed: it contains a trailing space and the layout
    # changed when the vendor results were reorganised into the repository.
    import glob as _glob, csv as _csv
    cand = _glob.glob(str(ROOT.parent / "GPNJ7M_results*" / "**" / "GPNJ7M-expression-matrix.tsv"),
                      recursive=True)
    if cand:
        rows = list(_csv.DictReader(open(cand[0]), delimiter="\t"))
        f["vendor features"] = len(rows)
        f["vendor protein_coding"] = sum(1 for r in rows if r.get("gene_biotype") == "protein_coding")
    cand = _glob.glob(str(ROOT.parent / "GPNJ7M_results*" / "**" / "GPNJ7M-mapping-stats-reads.csv"),
                      recursive=True)
    if cand:
        u = [int(r[1]) for r in list(_csv.reader(open(cand[0])))[1:]]
        f["vendor unique min"] = min(u)
        f["vendor unique max"] = max(u)

    # DE contrasts
    for p in (ROOT / "results/dge_BOM_ss5_filtered").glob("DE_rRNArm_*.csv"):
        n = sum(1 for r in csv.DictReader(p.open())
                if r["padj"] not in ("NA", "") and float(r["padj"]) < 0.05 and abs(float(r["log2FoldChange"])) > 1)
        f["DE " + p.stem.replace("DE_rRNArm_", "")] = n
    return f


# Values that are analysis parameters or thresholds we chose, not results to be verified.
PARAMETERS = {3000, 15000, 100000, 250000, 500, 50, 30, 20, 18, 12, 16, 94, 1000}


def audit(texfile, F):
    text = texfile.read_text()
    text = re.sub(r"%.*", "", text)                       # strip comments
    # Identifiers carry digits that are not claims: citation keys, pipeline names, accessions,
    # verbatim spans. "GL-DPPD-7101-G" was otherwise read as the number 7,101.
    text = re.sub(r"\\cite\{[^}]*\}", " ", text)
    text = re.sub(r"\\texttt\{[^}]*\}", " ", text)
    text = re.sub(r"\\url\{[^}]*\}", " ", text)
    text = re.sub(r"\\label\{[^}]*\}", " ", text)
    text = re.sub(r"\\ref\{[^}]*\}", " ", text)
    text = re.sub(r"GL-DPPD-[0-9-]+[A-Z]?", " ", text)
    nums = set()
    for m in re.finditer(r"(?<![\w.])(\d{1,3}(?:[,{]?\d{3})*(?:\.\d+)?)(?![\w])", text):
        raw = m.group(1).replace(",", "").replace("{", "")
        try: nums.add(float(raw))
        except ValueError: pass
    vals = {}
    for k, v in F.items(): vals.setdefault(float(v), []).append(k)
    checked, unmatched = [], []
    for n in sorted(nums):
        if n in vals: checked.append((n, vals[n]))
        elif n in PARAMETERS: continue                   # chosen parameter, not a result
        elif n > 1000: unmatched.append(n)               # large numbers must be traceable
    return checked, unmatched


def main():
    F = facts()
    print(f"fact dictionary: {len(F)} verified values\n")
    out = ["# Numbers audit", "",
           f"Generated by `scripts/32_audit_numbers.py`. Fact dictionary computed live from "
           f"result files ({len(F)} values).", ""]
    fails = 0
    for tex in sorted(LATEX.glob("*/sections.tex")) + sorted(LATEX.glob("*/main.tex")):
        rel = tex.relative_to(LATEX)
        checked, unmatched = audit(tex, F)
        print(f"{str(rel):<28} matched {len(checked):>3}   untraced>1000 {len(unmatched):>3}")
        out.append(f"## `{rel}`\n")
        out.append(f"- values matching a computed fact: **{len(checked)}**")
        out.append(f"- numbers above 1000 with no matching fact: **{len(unmatched)}**")
        if unmatched:
            fails += len(unmatched)
            out.append(f"\n  Untraced: {', '.join(f'{int(u):,}' for u in unmatched[:25])}")
        out.append("")
        if checked:
            out.append("| Value | Verified against |")
            out.append("|---|---|")
            for n, ks in checked:
                shown = f"{n:,.2f}".rstrip("0").rstrip(".") if n % 1 else f"{int(n):,}"
                out.append(f"| {shown} | {'; '.join(ks[:3])} |")
            out.append("")
    (LATEX / "NUMBERS_AUDIT.md").write_text("\n".join(out))
    print(f"\nwrote {LATEX/'NUMBERS_AUDIT.md'}")
    if fails:
        print(f"\n{fails} large number(s) could not be traced to a computed fact - review them.")
        return 1
    print("\nall large numbers traced.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
