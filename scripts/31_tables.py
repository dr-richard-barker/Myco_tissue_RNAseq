#!/usr/bin/env python3
"""Emit manuscript tables as booktabs LaTeX fragments, plus large supplementary data as CSV.

Every value is read from a result file; nothing is transcribed by hand. Fragments are
\\input{} from the manuscripts so the same numbers appear in all three documents and cannot
drift between them.
"""
import csv, pathlib, re, statistics, collections, sys

ROOT = pathlib.Path(__file__).resolve().parents[1]
OUT = ROOT / "latex" / "tables"; OUT.mkdir(parents=True, exist_ok=True)
esc = lambda s: str(s).replace("_", r"\_").replace("%", r"\%").replace("&", r"\&")
num = lambda n: f"{int(n):,}"


def write(name, body):
    (OUT / f"{name}.tex").write_text(body)
    print(f"  {name}.tex")


def tabular(cols, header, rows, caption, label, notes=None, small=True):
    b = []
    b.append(r"\begin{table}[htbp]")
    b.append(r"\caption{" + caption + r"}\label{" + label + r"}")
    if small: b.append(r"\footnotesize")
    b.append(r"\begin{tabular}{@{}" + cols + r"@{}}")
    b.append(r"\toprule")
    b.append(" & ".join(header) + r" \\")
    b.append(r"\midrule")
    for r in rows: b.append(" & ".join(str(x) for x in r) + r" \\")
    b.append(r"\bottomrule")
    b.append(r"\end{tabular}")
    if notes:
        b.append(r"\begin{tablenotes}[flushleft]\footnotesize\item " + notes + r"\end{tablenotes}")
    b.append(r"\end{table}")
    return "\n".join(b) + "\n"


# ---------------- T1: read budget ----------------
bud = list(csv.DictReader((ROOT / "results/read_budget_BOM_ss5.csv").open()))
order = ["Exuding mycelium", "Fuzzy mycelium", "Exudophore", "Nodule"]
bud.sort(key=lambda r: (order.index(r["tissue"]), r["well"]))
rows = []
for r in bud:
    keep = "" if int(r["mRNA"]) >= 15000 else r"$\dagger$"
    rows.append([esc(r["tissue"]), r["well"], num(r["raw"]), num(r["rRNA"]),
                 num(r["mRNA"]), f'{float(r["pct"]):.2f}', num(r["det10"]), keep])
write("T1_read_budget", tabular(
    "llrrrrrc",
    ["Tissue", "Well", "Raw reads", "rRNA", "mRNA", r"mRNA \%", r"Genes $\geq$10", ""],
    rows,
    r"Per-library read budget against the BOM\_ss5 reference after UMI deduplication. "
    r"rRNA and mRNA are counts assigned to features of the respective biotype.",
    "tab:budget",
    r"$\dagger$ excluded from the differential expression analysis "
    r"(fewer than 15,000 assigned non-rRNA counts)."))

# ---------------- T2: reference comparison ----------------
tm = collections.defaultdict(list)
for p in (ROOT / "qc/testmap").glob("*.summary.txt"):
    t = p.read_text()
    g = lambda rx: float(re.search(rx, t).group(1))
    tot = g(r"Total reads: *(\d+)")
    tm[p.name.split("__")[0]].append((100*g(r"Aligned 1 time: *(\d+)")/tot,
                                      100*g(r"Aligned >1 times: *(\d+)")/tot,
                                      100*g(r"Aligned 0 time: *(\d+)")/tot))
meta = {"PC9": ("GCF\\_014466165.1", "contig", "11,849"),
        "PC9.15": ("GCA\\_029852705.2", "chromosome", "13,556"),
        "BOM\\_ss5": ("GCA\\_056149245.1", "contig", "12,705"),
        "BOM\\_ss14": ("GCA\\_056149315.1", "contig", "13,310")}
rows = []
for k in sorted(tm, key=lambda x: -statistics.mean(v[0] for v in tm[x])):
    e = esc(k); u, m, un = (statistics.mean(v[i] for v in tm[k]) for i in range(3))
    acc, lvl, ng = meta.get(e, ("", "", ""))
    rows.append([e, acc, lvl, ng, f"{u:.2f}", f"{m:.2f}", f"{un:.2f}"])
write("T2_references", tabular(
    "lllrrrr",
    ["Assembly", "Accession", "Level", "Genes", r"Unique \%", r"Multi \%", r"Unaligned \%"],
    rows,
    r"Reference candidates compared empirically. 250{,}000 reads per library mapped to "
    r"genome-only indices; values are means across all 16 libraries.",
    "tab:refs",
    r"The RefSeq-designated reference genome for the species (PC9) performed worst."))

# ---------------- T3: FDR-surviving modules ----------------
conf = [r for r in csv.DictReader((ROOT / "results/wgcna/rrna_confound_test.csv").open())
        if float(r["fdr_tissue"]) < 0.05 and abs(float(r["r_tissue"])) > 0.7]
titles = {r["module"]: r for r in csv.DictReader((OUT.parent.parent / "results/wgcna/module_titles.tsv").open(), delimiter="\t")}
rows = []
for r in sorted(conf, key=lambda x: float(x["fdr_tissue"])):
    t = titles.get(r["module"], {})
    ttl = t.get("title", "")
    ttl = (ttl[:52] + "...") if len(ttl) > 55 else ttl
    rows.append([esc(r["module"]), t.get("size", ""), esc(r["tissue"]),
                 f'{float(r["r_tissue"]):+.2f}', f'{float(r["fdr_tissue"]):.4f}',
                 f'{float(r["r_rRNAfrac"]):+.2f}', f'{float(r["r_tissue_partial"]):+.2f}',
                 esc(ttl)])
write("T3_modules", tabular(
    "lrlrrrrp{42mm}",
    ["Module", "Genes", "Tissue", "$r$", "FDR", r"$r_{\mathrm{rRNA}}$", r"$r_{\mathrm{partial}}$",
     "Enriched terms"],
    rows,
    r"Co-expression modules whose tissue association survives Benjamini--Hochberg correction "
    r"across all 168 module\,$\times$\,tissue tests.",
    "tab:modules",
    r"$r_{\mathrm{rRNA}}$, correlation of the module eigengene with per-library rRNA "
    r"fraction; $r_{\mathrm{partial}}$, tissue correlation controlling for it."))

# ---------------- Supplementary ----------------
rs = list(csv.DictReader((ROOT / "metadata/runsheet.csv").open()))
rows = [[esc(r["sample_name"]), r["well"], esc(r["Factor Value[Tissue]"])] for r in rs]
write("S1_samples", tabular("lll", ["Library", "Well", "Tissue"], rows,
    r"Sample manifest. All libraries were sequenced on one flowcell and lane.", "tab:s1"))

de_rows = []
for f in sorted((ROOT / "results/dge_BOM_ss5_filtered").glob("DE_rRNArm_*.csv")):
    n = sum(1 for r in csv.DictReader(f.open())
            if r["padj"] not in ("NA", "") and float(r["padj"]) < 0.05
            and abs(float(r["log2FoldChange"])) > 1)
    a, b = f.stem.replace("DE_rRNArm_", "").split("_vs_")
    de_rows.append([esc(a.replace(".", " ")), esc(b.replace(".", " ")), n])
de_rows.sort(key=lambda r: -r[2])
write("S2_de", tabular("llr", ["Tissue A", "Tissue B", "DE genes"], de_rows,
    r"Differential expression across all six pairwise contrasts (12 retained libraries, "
    r"rRNA-removed track; $p_{\mathrm{adj}}<0.05$, $|\log_2\mathrm{FC}|>1$).", "tab:s2"))

gem = list(csv.DictReader((ROOT / "results/figure_data/gem_stats.csv").open()))
lbl = {"draft": "EC-mapped draft", "medium": "Transport + MNM medium", "gapfilled": "Targeted gapfilling"}
write("S7_gem", tabular("lrrrrr",
    ["Stage", "Reactions", "Metabolites", "Genes", "Carry flux", "Gapfilled"],
    [[lbl[g["stage"]], num(g["reactions"]), num(g["metabolites"]), num(g["genes"]),
      num(g["carrying"]), g["gapfilled"]] for g in gem],
    r"Genome-scale metabolic reconstruction for \textit{P.\ ostreatus} BOM\_ss5.", "tab:s7",
    r"Gapfilled reactions carry no gene association."))

# large tables exported as CSV supplementary data
import shutil
for src, dest in [("results/wgcna/module_titles.tsv", "S4_module_titles.tsv"),
                  ("results/wgcna/module_enrichment.tsv", "S5_module_enrichment.tsv"),
                  ("results/annotation/bom_ss5_cazymes.tsv", "S6_cazymes.tsv"),
                  ("results/tissue_models_BOM_ss5/markers_robust.csv", "S3_robust_markers.csv"),
                  ("results/annotation/bom_ss5_secretome.tsv", "S8_secretome.tsv")]:
    shutil.copy(ROOT / src, OUT / dest); print(f"  {dest} (data file)")

print(f"\nwrote {OUT}")
