#!/usr/bin/env python3
"""
export_web_assets.py
Extracts project results into JSON data files, converts PDF figures to PNG,
and stages manuscript PDFs into docs/ for GitHub Pages.
"""

import os
import sys
import glob
import json
import csv
import shutil
import subprocess

WORKSPACE = "/Users/drb_laptop/Documents/Myco_tissue_RNAseq"
DOCS_DIR = os.path.join(WORKSPACE, "docs")
DATA_DIR = os.path.join(DOCS_DIR, "data")
FIG_DIR = os.path.join(DOCS_DIR, "figures")
PAPER_DIR = os.path.join(DOCS_DIR, "papers")

os.makedirs(DATA_DIR, exist_ok=True)
os.makedirs(FIG_DIR, exist_ok=True)
os.makedirs(PAPER_DIR, exist_ok=True)

print("1. Converting LaTeX PDF figures to 250 DPI PNGs...")
pdf_figures = sorted(glob.glob(os.path.join(WORKSPACE, "latex/figures/*.pdf")))
print(f"Found {len(pdf_figures)} figure PDFs.")

for pdf_path in pdf_figures:
    base_name = os.path.splitext(os.path.basename(pdf_path))[0]
    out_png = os.path.join(FIG_DIR, f"{base_name}.png")
    cmd = [
        "gs",
        "-dSAFER",
        "-dBATCH",
        "-dNOPAUSE",
        "-sDEVICE=png16m",
        "-r250",
        f"-sOutputFile={out_png}",
        pdf_path
    ]
    subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)
    if os.path.exists(out_png):
        size_kb = os.path.getsize(out_png) / 1024
        print(f"  Converted {base_name}.png ({size_kb:.1f} KB)")
    else:
        print(f"  Error converting {pdf_path}")

print("\n2. Staging manuscript draft PDFs...")
manuscripts = glob.glob(os.path.join(WORKSPACE, "manuscript/*.pdf"))
for m in manuscripts:
    dest = os.path.join(PAPER_DIR, os.path.basename(m))
    shutil.copy2(m, dest)
    print(f"  Copied {os.path.basename(m)}")

print("\n3. Generating docs/data/read_budget.json...")
# Read budget from results/read_budget_BOM_ss5.csv and flag low from sample_yield.csv
flagged_samples = set()
yield_csv = os.path.join(WORKSPACE, "results/dge_BOM_ss5_filtered/sample_yield.csv")
if os.path.exists(yield_csv):
    with open(yield_csv, mode="r", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            if row.get("flagged_low") == "TRUE":
                flagged_samples.add(row["sample"])

budget_csv = os.path.join(WORKSPACE, "results/read_budget_BOM_ss5.csv")
read_budget_list = []
with open(budget_csv, mode="r", encoding="utf-8") as f:
    reader = csv.DictReader(f)
    for row in reader:
        sample_id = row["sample"]
        is_flagged = sample_id in flagged_samples
        read_budget_list.append({
            "sample": sample_id,
            "well": row["well"],
            "tissue": row["tissue"],
            "raw_reads": int(row["raw"]),
            "rrna_reads": int(row["rRNA"]),
            "mrna_reads": int(row["mRNA"]),
            "genes_det1": int(row["det1"]),
            "genes_det10": int(row["det10"]),
            "mrna_pct": round(float(row["pct"]), 2),
            "flagged_low": is_flagged,
            "status": "Excluded (<15k counts)" if is_flagged else "Retained"
        })

# Sort by well (1A..1H, 2A..2H)
def well_sort_key(item):
    w = item["well"]
    return (w[0], int(w[1:]))

read_budget_list.sort(key=lambda x: (x["well"][0], x["well"][1]))

with open(os.path.join(DATA_DIR, "read_budget.json"), "w", encoding="utf-8") as f:
    json.dump(read_budget_list, f, indent=2)
print("  Created read_budget.json")

print("\n4. Generating docs/data/pca_data.json...")
# All 16 PCA from results/dge_BOM_ss5/PCA_rRNArm.csv
pca_all_csv = os.path.join(WORKSPACE, "results/dge_BOM_ss5/PCA_rRNArm.csv")
pca_all = []
if os.path.exists(pca_all_csv):
    with open(pca_all_csv, mode="r", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            sample_id = row["sample"].strip('"')
            well = sample_id.split("_")[-1].replace("sample", "").strip()
            pca_all.append({
                "sample": sample_id,
                "well": well,
                "tissue": row["Tissue"].strip('"').replace(".", " "),
                "pc1": round(float(row["PC1"]), 2),
                "pc2": round(float(row["PC2"]), 2),
                "mrna": int(row["mrna"]),
                "flagged_low": sample_id in flagged_samples
            })

# 12 Filtered PCA from results/dge_BOM_ss5_filtered/PCA_rRNArm.csv
pca_filt_csv = os.path.join(WORKSPACE, "results/dge_BOM_ss5_filtered/PCA_rRNArm.csv")
pca_filtered = []
if os.path.exists(pca_filt_csv):
    with open(pca_filt_csv, mode="r", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            sample_id = row["sample"].strip('"')
            well = sample_id.split("_")[-1].replace("sample", "").strip()
            pca_filtered.append({
                "sample": sample_id,
                "well": well,
                "tissue": row["Tissue"].strip('"').replace(".", " "),
                "pc1": round(float(row["PC1"]), 2),
                "pc2": round(float(row["PC2"]), 2),
                "mrna": int(row["mrna"]),
                "flagged_low": False
            })

pca_payload = {
    "all_16": {
        "description": "All 16 libraries (unfiltered). PC1 strongly reflects library depth (cor = 0.681).",
        "cor_depth": 0.681,
        "points": pca_all
    },
    "filtered_12": {
        "description": "12 retained libraries after excluding 4 noise samples (<15k counts). PC1 is depth-independent (cor = -0.003).",
        "cor_depth": -0.003,
        "points": pca_filtered
    }
}
with open(os.path.join(DATA_DIR, "pca_data.json"), "w", encoding="utf-8") as f:
    json.dump(pca_payload, f, indent=2)
print("  Created pca_data.json")

print("\n5. Generating docs/data/robust_markers.json...")
markers_csv = os.path.join(WORKSPACE, "results/tissue_models_BOM_ss5/markers_robust.csv")
robust_markers = []
with open(markers_csv, mode="r", encoding="utf-8") as f:
    reader = csv.DictReader(f)
    for row in reader:
        robust_markers.append({
            "tissue": row["tissue"],
            "gene": row["gene"],
            "protein_name": row["protein_name"] if row["protein_name"] else "Uncharacterized protein",
            "tau": float(row["tau"])
        })

with open(os.path.join(DATA_DIR, "robust_markers.json"), "w", encoding="utf-8") as f:
    json.dump(robust_markers, f, indent=2)
print(f"  Created robust_markers.json ({len(robust_markers)} markers)")

print("\n6. Generating docs/data/modules.json...")
# WGCNA modules surviving Benjamini-Hochberg FDR correction
modules = [
    {
        "module": "greenyellow",
        "color": "#a3e635",
        "genes": 142,
        "tissue": "Exudophore",
        "r": -0.92,
        "fdr": 0.0034,
        "r_rrna": 0.28,
        "r_partial": -0.92,
        "enriched_terms": "Pexophagy / response to calcium ion / cellular homeostasis"
    },
    {
        "module": "lightcyan",
        "color": "#67e8f9",
        "genes": 99,
        "tissue": "Exudophore",
        "r": 0.91,
        "fdr": 0.0034,
        "r_rrna": -0.15,
        "r_partial": 0.91,
        "enriched_terms": "Maturation of SSU-rRNA / CENP-A containing chromatin organization"
    },
    {
        "module": "white",
        "color": "#e2e8f0",
        "genes": 55,
        "tissue": "Nodule",
        "r": -0.88,
        "fdr": 0.0091,
        "r_rrna": -0.10,
        "r_partial": -0.92,
        "enriched_terms": "One-carbon metabolic process / methyltransferase activity"
    },
    {
        "module": "lightsteelblue1",
        "color": "#93c5fd",
        "genes": 34,
        "tissue": "Exudophore",
        "r": 0.85,
        "fdr": 0.0155,
        "r_rrna": -0.36,
        "r_partial": 0.85,
        "enriched_terms": "Chromosome segregation / cell division"
    },
    {
        "module": "brown",
        "color": "#b45309",
        "genes": 320,
        "tissue": "Nodule",
        "r": 0.85,
        "fdr": 0.0155,
        "r_rrna": -0.09,
        "r_partial": 0.85,
        "enriched_terms": "Oxidoreductase activity / acting on paired donors with O2 incorporation"
    },
    {
        "module": "saddlebrown",
        "color": "#78350f",
        "genes": 49,
        "tissue": "Exudophore",
        "r": 0.84,
        "fdr": 0.0179,
        "r_rrna": 0.07,
        "r_partial": 0.88,
        "enriched_terms": "Methane metabolism / Microbial metabolism in diverse environments"
    }
]
with open(os.path.join(DATA_DIR, "modules.json"), "w", encoding="utf-8") as f:
    json.dump(modules, f, indent=2)
print("  Created modules.json")

print("\n7. Generating docs/data/figures_meta.json...")
figures_meta = [
    {
        "id": "fig1_read_fate",
        "number": "Figure 1",
        "title": "Read Fate & Ribosomal RNA Dominance",
        "category": "qc",
        "paper": "Paper 1 (Resource & QC)",
        "caption": "Fate of 92.4 million raw 3'-tag reads across 16 libraries. Nuclear rDNA accounts for 53.5% mean read share, while unannotated mitochondrial rRNA accounts for ~60% of reads surviving nuclear depletion. Only 2.7% of raw reads assign to protein-coding features."
    },
    {
        "id": "fig2_reference_rdna",
        "number": "Figure 2",
        "title": "Empirical Reference Selection & Nuclear rDNA Resolution",
        "category": "qc",
        "paper": "Paper 1 (Resource & QC)",
        "caption": "Comparison of candidate genome assemblies. BOM_ss5 uniquely resolves the tandem nuclear rDNA repeat array (JBQVBD010000012.1: 2.29-2.34 Mb) that PC9 and PC9.15 collapse into multi-mapping noise, yielding +6.9% higher assigned counts across all 16 samples."
    },
    {
        "id": "fig3_pca_retention",
        "number": "Figure 3",
        "title": "PCA & Elimination of the Depth Artefact",
        "category": "de",
        "paper": "Paper 1 (Resource & QC)",
        "caption": "Principal Component Analysis before and after filtering four noise-dominated libraries (<15,000 counts: 1C, 2E, 2B, 2A). Excluding these four reduces cor(PC1, depth) from +0.681 to -0.003 and triples the number of differentially expressed genes from 76 to 281."
    },
    {
        "id": "fig4_wgcna",
        "number": "Figure 4",
        "title": "WGCNA Module-Trait Associations & FDR Significance",
        "category": "de",
        "paper": "Paper 1 & Paper 2",
        "caption": "Co-expression network analysis. Only Nodule and Exudophore yield modules that survive Benjamini-Hochberg FDR correction (FDR < 0.05). Fuzzy and Exuding mycelium produce zero significant modules."
    },
    {
        "id": "fig5_exudophore_markers",
        "number": "Figure 5",
        "title": "Transcriptional Profile of the Exudophore & Nodule",
        "category": "tissue",
        "paper": "Paper 2 (Exudophore Biology)",
        "caption": "Expression profiles of top robust marker genes across tissues. The exudophore exhibits an oxidative and secretory profile centered on alcohol oxidase 1 (AOX1), glyoxal oxidase (GLOX), acetate--CoA ligase, and LPMO AA9. Nodule is dominated by cerato-platanin (CP1) and small HSPs."
    },
    {
        "id": "fig6_tissue_models",
        "number": "Figure 6",
        "title": "Genome-Scale Metabolic Model Architecture & Coverage",
        "category": "gem",
        "paper": "Paper 3 (Systems Metabolism)",
        "caption": "Overview of the Pleurotus ostreatus BOM_ss5 draft metabolic model, displaying protein annotation coverage, reaction distribution across metabolic subsystems, and reaction activity states."
    },
    {
        "id": "fig7_reaction_expression",
        "number": "Figure 7",
        "title": "Reaction-Level Flux & Expression Contrasts Across Tissues",
        "category": "gem",
        "paper": "Paper 3 (Systems Metabolism)",
        "caption": "Mapping tissue-specific gene expression to GPR (gene-protein-reaction) rules in the metabolic model, highlighting differential activation of one-carbon oxidation, polyol synthesis, and secondary metabolism in the exudophore."
    },
    {
        "id": "fig8_mitogenome_map",
        "number": "Figure 8",
        "title": "Complete Mitogenome Architecture of Pleurotus ostreatus",
        "category": "mito",
        "paper": "Paper 4 (Mitogenome)",
        "caption": "Annotation of the 73.2 kb mitochondrial contig CM148777.1, showing 13 core oxidative phosphorylation genes, 25 tRNAs, large and small subunit rRNAs, and multiple intron-encoded open reading frames."
    },
    {
        "id": "fig9_introns",
        "number": "Figure 9",
        "title": "Mitochondrial Group I/II Intron & Homing Endonuclease Distribution",
        "category": "mito",
        "paper": "Paper 4 (Mitogenome)",
        "caption": "Catalogue of mitochondrial introns and encoded LAGLIDADG / GIY-YIG homing endonuclease genes (HEGs) within cox1, cob, and rrnL genes, demonstrating their role in mitogenome expansion."
    },
    {
        "id": "fig10_mito_expression",
        "number": "Figure 10",
        "title": "Mitochondrial RNA-Seq Expression Dynamics",
        "category": "mito",
        "paper": "Paper 4 (Mitogenome)",
        "caption": "Coverage profiles and expression quantification across the mitotranscriptome, showing robust transcription of core respiratory subunits and co-transcription of intron-encoded ORFs."
    },
    {
        "id": "fig11_circular_map",
        "number": "Figure 11",
        "title": "Circular Mitogenome Map of Pleurotus ostreatus",
        "category": "mito",
        "paper": "Paper 4 & Paper 5",
        "caption": "High-resolution circular visualization of the complete P. ostreatus mitogenome, displaying strand-specific coding sequences, GC content, GC skew, and structural rRNA/tRNA coordinates."
    },
    {
        "id": "fig12_comparative_size",
        "number": "Figure 12",
        "title": "Comparative Mitogenomics Across Pleurotus & Agaricales",
        "category": "mito",
        "paper": "Paper 5 (Comparative Mitogenomics)",
        "caption": "Mitogenome size variation across twelve Pleurotus species and isolates (60.7 kb to 103.0 kb), showing a near-perfect correlation (r = +0.96) between genome size and open reading frame count."
    },
    {
        "id": "fig13_orfs_retrograde",
        "number": "Figure 13",
        "title": "Mitochondrial Duplicated ORFs & Nuclear Retrograde Signaling",
        "category": "mito",
        "paper": "Paper 5 (Comparative Mitogenomics)",
        "caption": "Identification of a duplicated uncharacterized mitochondrial ORF and expression dynamics of the 423 nuclear-encoded mitochondrial proteome genes indicating coordinated retrograde signaling."
    },
    {
        "id": "figS1_soft_threshold",
        "number": "Figure S1",
        "title": "WGCNA Soft-Threshold Power Diagnostics",
        "category": "supp",
        "paper": "Paper 1 & Paper 2",
        "caption": "Scale-free topology fit index and mean connectivity across soft-threshold powers beta, justifying the selection of beta = 14 for signed network construction."
    },
    {
        "id": "figS2_utr_extension",
        "number": "Figure S2",
        "title": "3' UTR Extension Benchmark",
        "category": "supp",
        "paper": "Paper 1 (Resource & QC)",
        "caption": "Benchmark evaluating read assignment gain from neighbour-capped 3' end extension. Extending 97.7% of gene models shifted assignment from 15.5% to 15.6% (+0.1%), proving truncated UTRs were not the bottleneck."
    },
    {
        "id": "figS3_tau_robustness",
        "number": "Figure S3",
        "title": "Tissue Specificity (Tau) Distribution & Replicate Robustness Filter",
        "category": "supp",
        "paper": "Paper 1 & Paper 2",
        "caption": "Distribution of tau tissue-specificity indices and the pass rate of top-50 candidate markers under the stringent across-replicate robustness test: Exudophore (74%), Nodule (44%), Fuzzy (2%), Exuding (2%)."
    },
    {
        "id": "figS4_cazy_secretome",
        "number": "Figure S4",
        "title": "CAZyme Classification & Predicted Secretome Repertoire",
        "category": "supp",
        "paper": "Paper 1 & Paper 3",
        "caption": "Distribution of 423 Carbohydrate-Active Enzymes (CAZymes) across GH, GT, PL, CE, AA, and CBM classes, alongside Swiss-Prot homology-transferred secretome annotations."
    },
    {
        "id": "figS5_gem_stages",
        "number": "Figure S5",
        "title": "Metabolic Reconstruction Progression Stages",
        "category": "supp",
        "paper": "Paper 3 (Systems Metabolism)",
        "caption": "Evolution of the metabolic model from raw draft (81.5% blocked) to Mycoponic Nutrient Medium (MNM v3) integration, and targeted backward MILP gapfilling (achieving 125.0 biomass flux)."
    },
    {
        "id": "figS6_duplication",
        "number": "Figure S6",
        "title": "UMI Deduplication Rates Across Libraries",
        "category": "supp",
        "paper": "Paper 1 (Resource & QC)",
        "caption": "Deduplication rates evaluated across 14 nt UMIs, demonstrating the vital necessity of directional UMI collapsing in 3'-tag single-end data."
    },
    {
        "id": "figS7_intergenic_peaks",
        "number": "Figure S7",
        "title": "Intergenic Coverage Peaks Pinpointing Hidden rDNA Loci",
        "category": "supp",
        "paper": "Paper 1 (Resource & QC)",
        "caption": "Clustering pooled unassigned reads revealed huge discrete coverage peaks corresponding to the unannotated tandem nuclear rDNA array and mitochondrial large subunit rRNA."
    }
]
with open(os.path.join(DATA_DIR, "figures_meta.json"), "w", encoding="utf-8") as f:
    json.dump(figures_meta, f, indent=2)
print("  Created figures_meta.json")

print("\nAsset generation complete!")
