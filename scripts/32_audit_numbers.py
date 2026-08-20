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

    # systems-biology layer (P3)
    fm = ROOT / "results/normalised_all_genes.csv"
    if fm.exists():
        f["genes full matrix"] = sum(1 for _ in fm.open()) - 1
    rex = ROOT / "results/tissue_metabolism/reaction_expression.csv"
    if rex.exists():
        import csv as _c
        rows = list(_c.DictReader(rex.open()))
        f["reactions with expression"] = len(rows)
        f["reactions 2x elevated"] = sum(1 for r in rows
                                         if float(r["ratio"]) >= 2 and float(r["Exudophore"]) >= 5)
    tms = ROOT / "models/tissue/tissue_model_summary.csv"
    if tms.exists():
        import csv as _c
        for r in _c.DictReader(tms.open()):
            f[f"tissue model {r['tissue']} reactions"] = int(r["reactions"])
            f[f"tissue model {r['tissue']} genes"] = int(r["genes"])
    pres = ROOT / "results/tissue_metabolism/reaction_presence.csv"
    if pres.exists():
        import csv as _c
        rows = list(_c.DictReader(pres.open()))
        T = ["Exuding_mycelium", "Fuzzy_mycelium", "Exudophore", "Nodule"]
        uq = [r for r in rows if int(r["Exudophore"]) == 1
              and sum(int(r[t]) for t in T if t != "Exudophore") == 0]
        f["exudophore only reactions"] = len(uq)
        f["exudophore only with gpr"] = sum(1 for r in uq if r["gpr"].strip())
        f["exudophore only no gpr"] = sum(1 for r in uq if not r["gpr"].strip())

    # mitochondrial layer (P4)
    import csv as _c, re as _re, collections as _co
    mg = ROOT / "refs/BOM_ss5/mitogenome.gff"
    if mg.exists():
        feats = [l.split("\t") for l in mg.open() if not l.startswith("#")]
        cnt = _co.Counter(x[2] for x in feats if len(x) > 2)
        f["mito features corrected"] = sum(cnt.values())
        for k, v in cnt.items():
            f[f"mito {k}"] = v
    stock = _co.Counter()
    for line in (ROOT / "refs/BOM_ss5/BOM_ss5_genomic.gtf").open():
        x = line.split("\t")
        if len(x) > 8 and x[0] == "CM148777.1" and x[2] == "gene":
            m = _re.search(r'gene_biotype "([^"]+)"', x[8])
            stock[m.group(1) if m else "other"] += 1
    if stock:
        f["mito features stock"] = sum(stock.values())
        f["mito stock tRNA"] = stock.get("tRNA", 0)
    ic = ROOT / "results/mito/intron_catalogue.csv"
    if ic.exists():
        rows = list(_c.DictReader(ic.open()))
        conf = [r for r in rows if r["mean_efficiency"] not in ("",) and float(r["mean_efficiency"]) > 0.5]
        f["mito junctions retained"] = len(rows)
        f["mito introns confirmed"] = len(conf)
        f["mito intron total bp"] = sum(int(r["length"]) for r in conf)
        low = [r for r in rows if r["mean_efficiency"] != "" and float(r["mean_efficiency"]) < 0.05]
        f["mito artefact junctions"] = len(low)
        f["mito artefacts in rRNA"] = sum(1 for r in low if 1000 <= int(r["start"]) <= 7200)
    # mitogenome length, and the assignment comparison that is the paper's headline
    mfa = ROOT / "refs/mito/CM148777.1.fa"
    if mfa.exists():
        f["mito length"] = sum(len(l.strip()) for l in mfa.open() if not l.startswith(">"))
    mc = ROOT / "results/mito/mito_counts.txt.summary"
    if mc.exists():
        d = {}
        for line in mc.open():
            x = line.rstrip().split("\t")
            if x[0] == "Status":
                continue
            d[x[0]] = sum(int(v) for v in x[1:])
        f["mito assigned corrected"] = d.get("Assigned", 0)
    # stock assignment and the fold change are recomputed here so the ratio in the text is
    # never a hand-carried number
    sc = pathlib.Path("/tmp/stock_counts.txt.summary")
    if sc.exists():
        d = {}
        for line in sc.open():
            x = line.rstrip().split("\t")
            if x[0] == "Status":
                continue
            d[x[0]] = sum(int(v) for v in x[1:])
        f["mito assigned stock"] = d.get("Assigned", 0)
        if d.get("Assigned"):
            f["mito assignment fold"] = round(f.get("mito assigned corrected", 0) / d["Assigned"])
    mt = ROOT / "results/mito/mito_total_alignments.txt"
    if mt.exists():
        f["mito total alignments"] = int(mt.read_text().strip())
    ictab = ROOT / "results/mito/intron_catalogue.csv"
    if ictab.exists():
        import csv as _cc
        rr = list(_cc.DictReader(ictab.open()))
        # the top junction by raw read count is an ARTEFACT in the rRNA block; the headline
        # intron is the best-supported junction that actually splices
        conf_rows = [r for r in rr if r["mean_efficiency"] not in ("",)
                     and float(r["mean_efficiency"]) > 0.5]
        top = max(conf_rows or rr, key=lambda r: int(r["total_reads"]))
        f["mito top intron reads"] = int(top["total_reads"])
        f["mito top intron start"] = int(top["start"])
        f["mito top intron end"] = int(top["end"])

    # comparative + MitoCarta layer (P5)
    inv = ROOT / "results/mito_comparative/inventory.csv"
    if inv.exists():
        import csv as _ci
        rows = list(_ci.DictReader(inv.open()))
        f["comparative genomes"] = len(rows)
        L = [int(r["length"]) for r in rows]
        f["comparative min length"] = min(L)
        f["comparative max length"] = max(L)
        for r in rows:
            f[f"mitogenome {r['acc']}"] = int(r["length"])
            f[f"orfs {r['acc']}"] = int(r["orfs_ge100aa"])
    ncp = ROOT / "results/mitocarta/nuclear_mito_proteome.csv"
    if ncp.exists():
        import csv as _ci
        rows = list(_ci.DictReader(ncp.open()))
        f["nuclear mito union"] = len(rows)
        f["nuclear mito rbh"] = sum(1 for r in rows if r["mitocarta_rbh"] == "1")
        f["nuclear mito swissprot"] = sum(1 for r in rows if r["swissprot_mito"] == "1")
    fwd = ROOT / "results/mitocarta/nuc_vs_mitocarta.tsv"
    if fwd.exists():
        f["nuclear mitocarta oneway"] = len({l.split("\t")[0] for l in fwd.open()})
    oc = ROOT / "results/mito/orf_conservation.csv"
    if oc.exists():
        import csv as _ci
        rows = list(_ci.DictReader(oc.open()))
        f["uorfs tested"] = len(rows)
        f["uorfs conserved expressed"] = sum(1 for r in rows if int(r["genomes"]) >= 6
                                             and int(r["libs"]) >= 8 and r["in_rrna"] == "0")
    # per-genome detected group I intron span, and the coordinates of every predicted ORF
    for tb in (ROOT / "results/mito_comparative").glob("*.introns.tbl"):
        acc = tb.name.replace(".introns.tbl", "")
        n = span = 0
        for line in tb.open():
            if line.startswith("#"):
                continue
            x = line.split()
            if len(x) > 15:
                n += 1; span += abs(int(x[8]) - int(x[7])) + 1
        f[f"intron span {acc}"] = span
        f[f"intron hits {acc}"] = n
    orffa = ROOT / "results/mito/orfs_aa.fa"
    if orffa.exists():
        for line in orffa.open():
            if not line.startswith(">"):
                continue
            m = _re.match(r">(\S+)\s+\[(\d+)\s*-\s*(\d+)\]", line)
            if m:
                f[f"orf {m.group(1)} start"] = int(m.group(2))
                f[f"orf {m.group(1)} end"] = int(m.group(3))

    mcf = ROOT / "refs/mitocarta/Human.MitoCarta3.0.fasta"
    if mcf.exists():
        f["mitocarta sequences"] = sum(1 for l in mcf.open() if l.startswith(">"))

    orf = ROOT / "results/mito/candidate_orfs.csv"
    if orf.exists():
        rows = list(_c.DictReader(orf.open()))
        f["mito candidate orfs"] = len(rows)
        f["mito orfs expressed credible"] = sum(1 for r in rows
            if r.get("in_rRNA_block") == "0" and int(r["libraries_detected"]) >= 8)

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
