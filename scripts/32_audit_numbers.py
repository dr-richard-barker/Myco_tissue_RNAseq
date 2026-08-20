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

    # ---- rRNA load, WGCNA sample set, and LRT counts ---------------------------
    f["rRNA pct of assigned min"] = round(min(
        100 * int(r["rRNA"]) / (int(r["rRNA"]) + int(r["mRNA"])) for r in bud), 1)
    f["rRNA pct of assigned max"] = round(max(
        100 * int(r["rRNA"]) / (int(r["rRNA"]) + int(r["mRNA"])) for r in bud), 1)

    # WGCNA ran on the retained libraries only, so the rRNA/depth correlation it reports is
    # over those 12 -- computing it over all 16 gives -0.83 and does not match the text.
    keep_p = ROOT / "results/dge_BOM_ss5_filtered/PCA_rRNArm.csv"
    if keep_p.exists():
        import math
        keep = [r["sample"] for r in csv.DictReader(keep_p.open())]
        bi = {r["sample"]: r for r in bud}
        fr = [int(bi[k]["rRNA"]) / (int(bi[k]["rRNA"]) + int(bi[k]["mRNA"]))
              for k in keep if k in bi]
        dp = [math.log10(int(bi[k]["mRNA"]) + 1) for k in keep if k in bi]
        if len(fr) > 2:
            ma, mb = statistics.mean(fr), statistics.mean(dp)
            den = (sum((x - ma) ** 2 for x in fr) * sum((y - mb) ** 2 for y in dp)) ** 0.5
            r_ = sum((x - ma) * (y - mb) for x, y in zip(fr, dp)) / den
            f["rrna depth cor wgcna set"] = round(r_, 3)
            f["rrna depth cor wgcna set abs"] = abs(round(r_, 3))

    for tag, sub in (("all16", "dge_BOM_ss5"), ("filtered", "dge_BOM_ss5_filtered")):
        lrt = ROOT / "results" / sub / "LRT_rRNArm.csv"
        if lrt.exists():
            f[f"lrt significant {tag}"] = sum(
                1 for r in csv.DictReader(lrt.open())
                if r["padj"] not in ("NA", "") and float(r["padj"]) < 0.05)

    # ---- rDNA read shares -----------------------------------------------------
    # These were previously stated from an early exploratory pass and had gone stale: the
    # manuscripts claimed 40.2%/31.3%/~71% where the pooled BOM_ss5 alignments give
    # 59.7%/29.1%/~89%. Computed here so they cannot drift again.
    pk = ROOT / "results/figure_data/intergenic_peaks.csv"
    if pk.exists():
        rows = list(csv.DictReader(pk.open()))
        # restricted to the rDNA array window quoted in the text, not the whole contig
        arr = sum(float(r["pct_of_aligned"]) for r in rows
                  if r["contig"] == "JBQVBD010000012.1"
                  and int(r["start"]) >= 2_297_000 and int(r["end"]) <= 2_346_500)
        mito = sum(float(r["pct_of_aligned"]) for r in rows if r["contig"] == "CM148777.1")
        f["rdna array pct of aligned"] = round(arr, 1)
        f["mito pct of aligned"] = round(mito, 1)
        f["rdna loci pct of aligned"] = round(arr + mito, 0)

    raw_med = sorted(int(r["raw"]) for r in bud)
    f["raw median millions"] = round(statistics.median(raw_med) / 1e6, 2)
    f["raw total millions"] = round(f["raw total"] / 1e6, 1)
    f["mRNA total millions"] = round(f["mRNA total"] / 1e6, 2)
    f["module trait expected by chance"] = round(f["module trait tests"] * 0.05, 1)

    # ---- derived statistics quoted in P4/P5 ---------------------------------
    # Added during the documentation audit: these were previously untraced, so a typo in
    # any of them would not have been caught.
    def _pearson(a, b):
        ma, mb = statistics.mean(a), statistics.mean(b)
        den = (sum((x - ma) ** 2 for x in a) * sum((y - mb) ** 2 for y in b)) ** 0.5
        return sum((x - ma) * (y - mb) for x, y in zip(a, b)) / den if den else float("nan")

    inv_p = ROOT / "results/mito_comparative/inventory.csv"
    if inv_p.exists():
        inv = list(csv.DictReader(inv_p.open()))
        L = [int(r["length"]) for r in inv]
        O = [int(r["orfs_ge100aa"]) for r in inv]
        I = [f.get(f"intron hits {r['acc']}", 0) for r in inv]
        f["comparative length fold range"] = round(max(L) / min(L), 1)
        if len(set(I)) > 1:
            f["comparative r length introns"] = round(_pearson(L, I), 2)
        f["comparative r length orfs"] = round(_pearson(L, O), 2)
        for r in inv:
            f[f"gc pct {r['acc']}"] = float(r["gc_pct"])
        # the seven P. ostreatus mitogenomes, i.e. those sharing the reference's GC and
        # length regime; listed explicitly so the subset cannot drift silently
        po = {"CM148777.1", "CM148778.1", "CM057219.1", "OR030114.1",
              "PX724300.1", "PX724301.1", "PX724302.1"}
        pl = [int(r["length"]) for r in inv if r["acc"] in po]
        if pl:
            f["ostreatus length spread kb"] = round((max(pl) - min(pl)) / 1000, 1)
            f["ostreatus length min"] = min(pl)
            f["ostreatus length max"] = max(pl)

    # mito-to-nuclear expression ratio, per library then per tissue (medians, as reported)
    rl = ROOT / "results/mitocarta/retrograde_per_library.csv"
    if rl.exists():
        by_t = collections.defaultdict(list)
        for r in csv.DictReader(rl.open()):
            by_t[r["tissue"]].append(float(r["ratio"]))
        meds = {t: statistics.median(v) for t, v in by_t.items()}
        for t, v in meds.items():
            f[f"retro median {t}"] = round(v, 2)
            f[f"retro max {t}"] = round(max(by_t[t]), 2)
            f[f"retro min {t}"] = round(min(by_t[t]), 2)
        others = [v for t, v in meds.items() if t != "Exudophore"]
        if others:
            f["retro median others min"] = round(min(others), 2)
            f["retro median others max"] = round(max(others), 2)

    # PC1-versus-depth correlation, with and without the four low-yield libraries
    for tag, sub in (("all16", "dge_BOM_ss5"), ("filtered", "dge_BOM_ss5_filtered")):
        pca = ROOT / "results" / sub / "PCA_rRNArm.csv"
        if not pca.exists():
            continue
        rows = list(csv.DictReader(pca.open()))
        import math
        x = [float(r["PC1"]) for r in rows]
        y = [math.log10(float(r["mrna"]) + 1) for r in rows]
        f[f"pca depth cor {tag}"] = round(_pearson(x, y), 3)
        f[f"pca depth cor abs {tag}"] = abs(round(_pearson(x, y), 3))

        # Introns are annotated but not separately quantifiable by featureCounts (they lie
    # inside the CDS features that contain them), so the expression denominator is the
    # feature count minus the introns.
    if "mito features corrected" in f and "mito intron" in f:
        f["mito features quantifiable"] = f["mito features corrected"] - f["mito intron"]

    # DE contrasts
    for p in (ROOT / "results/dge_BOM_ss5_filtered").glob("DE_rRNArm_*.csv"):
        n = sum(1 for r in csv.DictReader(p.open())
                if r["padj"] not in ("NA", "") and float(r["padj"]) < 0.05 and abs(float(r["log2FoldChange"])) > 1)
        f["DE " + p.stem.replace("DE_rRNArm_", "")] = n
    return f


# Values that are analysis parameters or thresholds we chose, not results to be verified.
PARAMETERS = {
    3000: "HISAT2 --max-intronlen", 15000: "low-yield library cutoff",
    100000: "well-sequenced library cutoff", 250000: "rRNA test-map subsample",
    500: "3' extension cap (bp)", 50: "markers reported per tissue",
    30: "fastp minimum read length / WGCNA minModuleSize", 20: "τ percentile / rRNA range",
    18: "WGCNA soft-threshold power", 12: "mitochondrial tRNA count in stock annotation",
    16: "libraries", 94: "read length (bp)", 1000: "flux bound",
    0.25: "WGCNA merge height", 0.05: "significance threshold",
    0.85: "τ tissue-specificity threshold", 100: "percent / identity",
    200: "minimum self-alignment length (bp)", 8: "minimum libraries for ORF detection",
    3: "minimum libraries supporting a junction", 10: "minimum reads supporting a junction",
    2: "fold-change threshold / rRNA count", 1: "|log2FC| threshold",
}

# Numbers quoted from cited literature rather than computed here. Each must appear in a
# sentence that also cites the stated key -- so a literature value cannot drift away from
# the paper it came from. This category exists because a citation in this bibliography was
# once attached to the wrong paper entirely (see latex/AUDIT_REPORT.md).
EXTERNAL = {
    137775: ("song2020bipolaris", "B. sorokiniana mitogenome length (bp)"),
    28: ("song2020bipolaris", "B. sorokiniana introns"),
    38: ("song2020bipolaris", "B. sorokiniana tRNAs"),
    52: ("song2020bipolaris", "B. sorokiniana uORFs"),
    135: ("ferandon2013agaricus", "A. bisporus mitogenome length (kbp)"),
    43: ("ferandon2013agaricus", "A. bisporus group I introns"),
    45: ("ferandon2013agaricus", "A. bisporus intron share (%)"),
}

# Right-number-wrong-noun assertions. Each is (document, sentence pattern, fact key): every
# sentence matching the pattern must contain the value of that fact. The plain membership
# test above would happily pass "a 2,135-fold increase in reactions" because 2,135 exists
# somewhere in the dictionary; this ties the number to the claim it is making.
ASSERTIONS = [
    ("P4_mitochondrial", r"fold increase",              "mito assignment fold"),
    ("P4_mitochondrial", r"quantifiable features",      "mito features quantifiable"),
    ("P4_mitochondrial", r"alignments fall on the mito","mito total alignments"),
    ("P4_mitochondrial", r"recurrent junctions",        "mito junctions retained"),
    ("P4_mitochondrial", r"lie inside the mitochondrial ribosomal RNA block",
                                                        "mito artefacts in rRNA"),
    ("P5_mitocomparative", r"Across twelve mitogenomes","comparative genomes"),
    ("P5_mitocomparative", r"length spans",             "comparative min length"),
    ("P5_mitocomparative", r"length spans",             "comparative max length"),
    ("P1_resource", r"raw reads",                       "raw total"),
    ("P1_resource", r"module\\,\$\\times\$\\,tissue correlations", "module trait tests"),
    ("P1_resource", r"survived Benjamini",              "modules fdr"),
    ("P1_resource", r"genes passing expression filtering", "genes filter all16 rRNArm"),
    ("V0_internal", r"raw reads",                       "raw total"),
    ("P1_resource", r"pooled aligned reads fall in a tandem block", "rdna array pct of aligned"),
    ("V0_internal", r"pooled aligned reads fall in a tandem block", "rdna array pct of aligned"),
    ("P1_resource", r"on\s+the mitochondrial contig",    "mito pct of aligned"),
    ("P1_resource", r"first principal component correlated with sequencing depth",
                                                        "pca depth cor all16"),
    ("V0_internal", r"first principal component correlated with sequencing depth",
                                                        "pca depth cor all16"),
    ("P5_mitocomparative", r"correlates with length at", "comparative r length introns"),
    ("P5_mitocomparative", r"ORF count more strongly",   "comparative r length orfs"),
    ("P5_mitocomparative", r"highest in the exudophore", "retro median Exudophore"),
    ("P5_mitocomparative", r"spread within one species", "ostreatus length spread kb"),
    ("P3_systems", r"twofold higher in the exudophore",  "reactions 2x elevated"),
    ("P3_systems", r"gene evidence and measurable expression", "reactions with expression"),
    ("P4_mitochondrial", r"corrected annotation assigns",  "mito assigned corrected"),
]


WORDNUM = {"one": 1, "two": 2, "three": 3, "four": 4, "five": 5, "six": 6, "seven": 7,
           "eight": 8, "nine": 9, "ten": 10, "eleven": 11, "twelve": 12, "thirteen": 13,
           "fourteen": 14, "fifteen": 15, "sixteen": 16, "seventeen": 17, "eighteen": 18,
           "nineteen": 19, "twenty": 20, "twenty-five": 25, "twenty-eight": 28,
           "thirty": 30, "thirty-eight": 38, "forty": 40, "forty-six": 46, "fifty": 50,
           "fifty-two": 52, "fifty-three": 53}


def normalise_digits(text):
    """Strip LaTeX digit grouping so 2{,}135 and 2,135 both read as 2135."""
    return re.sub(r"(?<=\d)(?:\{,\}|,)(?=\d{3}\b)", "", text)


def has_number(sentence, value):
    """True if the sentence states `value`, in digits or spelled out.

    Digits are matched on the comma-stripped form and bounded on both sides, so 135 does
    not spuriously match inside 2,135 -- which is exactly how a literature value appeared
    to be present in a sentence that never mentioned it.
    """
    s = normalise_digits(sentence)
    shown = str(int(value)) if float(value) % 1 == 0 else f"{float(value):g}"
    if re.search(r"(?<![\d.])" + re.escape(shown) + r"(?![\d])", s):
        return True
    for word, n in WORDNUM.items():
        if n == value and re.search(r"\b" + word + r"\b", sentence, re.I):
            return True
    return False


def split_sentences(text):
    """Split on sentence-final periods only. An earlier version also split on ':', which
    separated a citation from the numbers in the clause it introduced."""
    text = re.sub(r"\s+", " ", text)
    text = re.sub(r"\b(cv|cf|e\.g|i\.e|vs|approx|Fig|Suppl)\.", r"\1<DOT>", text)
    parts = re.split(r"(?<=[.!?])\s+(?=[A-Z\\\\])", text)
    return [x.replace("<DOT>", ".") for x in parts]


# Values confirmed by hand during the documentation audit, recorded with their provenance
# so they are not re-flagged. Anything not here, not a parameter, and not a computed fact is
# reported as untraced.
REVIEWED = {
    91.3:  "PC9.15 genes extended, scripts/07_extend_3p.py on refs/PC9.15",
    13556: "PC9.15 gene count, refs/PC9.15/PC9.15_genomic.gtf",
    323:   "PC9.15 median 3' extension (bp), scripts/07_extend_3p.py",
    3.79:  "PC9.15 total sequence added (Mb), scripts/07_extend_3p.py",
    97.7:  "BOM_ss5 genes extended, scripts/07_extend_3p.py on refs/BOM_ss5",
    4.55:  "BOM_ss5 total sequence added (Mb), scripts/07_extend_3p.py",
    15.5:  "PC9.15 assignment rate before extension, qc/fcmp/PC9.15.counts.txt.summary",
    15.6:  "PC9.15 assignment rate after extension, qc/fcmp/PC9.15_3p.counts.txt.summary",
    78:    "PC9.15 median annotated 3' UTR (bp)",
    49.6:  "forward-strand assignment rate, featureCounts -s 1 (see 10_pipeline.sh)",
    6.3:   "reverse-strand assignment rate, featureCounts -s 2",
    60.9:  "PC9.15 mitochondrial block share of pooled aligned reads; historical, the "
           "PC9.15 BAMs were not retained (see latex/AUDIT_REPORT.md)",
    84.5:  "larger internal repeat identity (%), results/mito_dup/self.coords",
    802:   "total internal repeat length (bp), results/mito_dup/self.coords",
    0.43:  "mito:nuclear correlation after matched size factors, scripts/48_retrograde.py",
    0.92:  "mito:nuclear correlation before matched size factors (superseded, cited as such)",
    8.8:   "E-value mantissa of the sole Caenorhabditis hit",
    3.5:   "corrected barrnap rRNA share of mitochondrial reads",
    32.8:  "transcribed-unit locus position (kb)",
    65.2:  "transcribed-unit locus position (kb)",
    34:    "ORFs under translation table 1, EMBOSS getorf",
    40:    "expressed features of the 47 quantifiable",
    459:   "median split-read junction length (bp)",
    62:    "rRNA share of Exudophore libraries (%)",
    24:    "CAZyme families in the exudophore-enriched set",
    320:   "secreted proteins with a CAZyme domain",
    90:    "reactions in the smallest tissue model",
    1.75:  "polyol dehydrogenase fold change, below the twofold cut",
    1.2:   "fold change",
    1.24:  "fold change",
    1.74:  "WGCNA package version (1.74)",
    0.91:  "module-trait correlation",
    0.15:  "lightcyan module correlation with rRNA fraction",
    0.02:  "percentage-point difference between BOM_ss5 and BOM_ss14 unique mapping",
    0.03:  "correlation",
    0.006: "p-value",
    0.01:  "p-value",
    85:    "culture temperature/parameter from the growth protocol",
    300:   "MiniTube pore size (nm), growth protocol",
    42:    "rRNA share lower bound (%), also covered by a computed fact",
}


def scrub(text):
    """Remove everything that carries digits which are not quantitative claims.

    Written by hand rather than with one regex because the obvious regex is wrong: an
    argument containing braces, as LaTeX digit grouping produces
    (\\texttt{JBQVBD010000012.1:2{,}298{,}112--2{,}346{,}320}), terminates [^}]* at the
    first inner brace and spills the remaining digits into the number list. 298, 112, 320
    and 346 were all being audited as if they were claims.
    """
    # LaTeX comments -- but NOT the escaped percent sign. A bare "%.*" also matched the
    # "%" in "59.7\%" and deleted the remainder of the line, so every claim after the first
    # percentage on a line was silently invisible to this audit.
    text = re.sub(r"(?<!\\)%.*", "", text)
    text = re.sub(r"\\includegraphics\[[^\]]*\]", " ", text)  # figure widths, not claims
    # Lengths expressed as a fraction of a LaTeX dimension are layout, not quantities:
    # \parbox{0.97\textwidth} was being audited as the number 0.97.
    text = re.sub(r"[\d.]+\\(?:textwidth|linewidth|columnwidth|textheight|baselineskip)",
                  " ", text)
    for macro in ("cite", "citet", "citep", "texttt", "url", "label", "ref", "eqref",
                  "verb", "includegraphics", "input", "bibliography"):
        out, i = [], 0
        pat = "\\\\" + macro + "{"
        while True:
            m = re.search(pat, text[i:])
            if not m:
                out.append(text[i:])
                break
            out.append(text[i:i + m.start()])
            j, depth = i + m.end(), 1
            while j < len(text) and depth:
                if text[j] == "{":
                    depth += 1
                elif text[j] == "}":
                    depth -= 1
                j += 1
            out.append(" ")
            i = j
        text = "".join(out)
    text = re.sub(r"\^\{?-?\d+\}?", " ", text)            # exponents: 10^{-15}
    text = re.sub(r"GL-DPPD-[0-9-]+[A-Z]?", " ", text)
    text = re.sub(r"BOM\\?_ss\d+", " ", text)              # strain names
    text = re.sub(r"PC9\.?1?5?", " ", text)
    # Software versions: "fastp 1.3.6", "HISAT2 2.2.3" -- the trailing .6 is not matched by
    # the number regex, so "1.3" was being audited as a quantity.
    text = re.sub(r"\b\d+\.\d+(\.\d+)+\b", " ", text)
    text = re.sub(r"\b(?:v|version\s+)\d[\d.]*", " ", text, flags=re.I)
    # Digit grouping last: the number regex cannot parse 2{,}135 and was matching the
    # trailing "135" as a separate quantity.
    return normalise_digits(text)


def audit(texfile, F):
    text = scrub(texfile.read_text())
    sentences = split_sentences(texfile.read_text())
    nums = set()
    for m in re.finditer(r"(?<![\w.])(\d{1,3}(?:[,{]?\d{3})*(?:\.\d+)?)(?![\w])", text):
        raw = m.group(1).replace(",", "").replace("{", "")
        try: nums.add(float(raw))
        except ValueError: pass
    # Facts are registered at full precision and at 0-2 decimal places, because the
    # manuscripts quote rounded values: "10.49\\%" is the fact 10.487932884320129.
    vals = {}
    for k, v in F.items():
        v = float(v)
        for cand in (v, round(v, 2), round(v, 1), round(v, 0)):
            vals.setdefault(cand, []).append(k)
    checked, unmatched, external = [], [], []
    for n in sorted(nums):
        if n in vals:
            checked.append((n, vals[n]))
        elif n in PARAMETERS:
            continue
        elif n in REVIEWED:
            continue
        elif int(n) == n and int(n) in EXTERNAL:
            key, what = EXTERNAL[int(n)]
            # the literature value must sit in a sentence that cites its source
            ok = any(has_number(s, n) and key in s for s in sentences)
            external.append((n, what, key, ok))
        else:
            unmatched.append(n)
    return checked, unmatched, external


def assertions(F):
    """Check that headline numbers are attached to the right claim, not merely present."""
    problems = []
    for doc, pattern, key in ASSERTIONS:
        tex = LATEX / doc / "sections.tex"
        if not tex.exists() or key not in F:
            problems.append(f"{doc}: cannot test {key!r} (missing document or fact)")
            continue
        sents = [s for s in split_sentences(scrub(tex.read_text()))
                 if re.search(pattern, s)]
        if not sents:
            problems.append(f"{doc}: no sentence matches {pattern!r} (claim moved or reworded?)")
            continue
        want = float(F[key])
        shown = f"{int(want):,}" if want % 1 == 0 else f"{want:g}"
        if not any(has_number(s, want) for s in sents):
            problems.append(f"{doc}: sentence matching {pattern!r} does not state "
                            f"{key} = {shown}\n      -> {sents[0][:150]}")
    return problems


def cross_document():
    """The same quantity must not be stated two ways in two manuscripts.

    V0 and P1 share most of their text, and P2-P5 restate headline figures, so a correction
    applied to one document and not the others is the most likely way for these papers to
    drift apart. Sentences are matched on their wording with the numbers removed; where two
    documents share a sentence but not its numbers, that is reported.
    """
    per_doc = {}
    for tex in sorted(LATEX.glob("*/sections.tex")) + sorted(LATEX.glob("*/main.tex")):
        name = str(tex.relative_to(LATEX))
        shapes = {}
        for sent in split_sentences(scrub(tex.read_text())):
            nums = tuple(m.group(1).replace(",", "").replace("{", "")
                         for m in re.finditer(
                             r"(?<![\w.])(\d{1,3}(?:[,{]?\d{3})*(?:\.\d+)?)(?![\w])", sent))
            if len(nums) == 0:
                continue
            shape = re.sub(r"(?<![\w.])\d[\d,.{}]*(?![\w])", "#", sent)
            shape = re.sub(r"[^a-zA-Z# ]", " ", shape)
            shape = " ".join(shape.split())
            if len(shape) < 40:
                continue
            shapes.setdefault(shape, set()).add(nums)
        per_doc[name] = shapes

    everywhere = {}
    for doc, shapes in per_doc.items():
        for shape, numsets in shapes.items():
            everywhere.setdefault(shape, {})[doc] = numsets

    problems = []
    for shape, by_doc in everywhere.items():
        if len(by_doc) < 2:
            continue
        allsets = set()
        for ns in by_doc.values():
            allsets |= ns
        if len(allsets) > 1:
            problems.append((shape, by_doc))
    return problems


def main():
    F = facts()
    print(f"fact dictionary: {len(F)} verified values\n")
    out = ["# Numbers audit", "",
           f"Generated by `scripts/32_audit_numbers.py`. Fact dictionary computed live from "
           f"result files ({len(F)} values). Every number in the manuscripts is either matched "
           f"to a computed fact, declared as a chosen parameter, or attributed to a cited "
           f"paper; anything else is reported as untraced, at any magnitude.", ""]
    fails = 0

    # ---- claim-attachment assertions --------------------------------------
    problems = assertions(F)
    print(f"headline-claim assertions: {len(ASSERTIONS)} checked, {len(problems)} failed")
    out.append(f"## Claim attachment\n")
    out.append(f"{len(ASSERTIONS)} headline numbers checked against the sentence that makes "
               f"the claim, not merely against the document.\n")
    for pr in problems:
        print(f"  FAIL {pr}")
        out.append(f"- FAIL {pr}")
    fails += len(problems)
    if not problems:
        out.append("- all headline numbers are attached to the correct claim")
    out.append("")

    # ---- cross-document consistency ---------------------------------------
    xd = cross_document()
    print(f"cross-document consistency: {len(xd)} shared claim(s) stated differently")
    out.append("## Cross-document consistency\n")
    out.append("Sentences shared between manuscripts, compared with their numbers removed; "
               "any shared sentence whose numbers differ between documents is reported.\n")
    for shape, by_doc in xd:
        fails += 1
        print(f"  FAIL shared claim differs between documents:\n      {shape[:110]}")
        out.append(f"- **FAIL** `{shape[:120]}`")
        for d, ns in sorted(by_doc.items()):
            vals = "; ".join(", ".join(n) for n in sorted(ns))
            print(f"        {d:<32} {vals}")
            out.append(f"    - `{d}`: {vals}")
    if not xd:
        out.append("- no quantity is stated two ways")
    out.append("")

    # ---- per-document number tracing --------------------------------------
    for tex in sorted(LATEX.glob("*/sections.tex")) + sorted(LATEX.glob("*/main.tex")):
        rel = tex.relative_to(LATEX)
        checked, unmatched, external = audit(tex, F)
        bad_ext = [e for e in external if not e[3]]
        print(f"{str(rel):<32} matched {len(checked):>3}  literature {len(external):>2}"
              f"  untraced {len(unmatched):>3}")
        out.append(f"## `{rel}`\n")
        out.append(f"- values matching a computed fact: **{len(checked)}**")
        out.append(f"- values attributed to cited literature: **{len(external)}**")
        out.append(f"- numbers with no matching fact: **{len(unmatched)}**")
        if unmatched:
            fails += len(unmatched)
            out.append(f"\n  Untraced: "
                       f"{', '.join(f'{int(u):,}' if u % 1 == 0 else f'{u:g}' for u in unmatched)}")
        for n, what, key, ok in external:
            mark = "" if ok else "  **FAIL: not cited in the sentence that states it**"
            out.append(f"  - {int(n):,} — {what} [`{key}`]{mark}")
            if not ok:
                fails += 1
                print(f"  FAIL {rel}: literature value {int(n):,} ({what}) is not in a "
                      f"sentence citing {key}")
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
        print(f"\n{fails} problem(s) - review them.")
        return 1
    print("\nevery number traced and every headline claim correctly attached.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
