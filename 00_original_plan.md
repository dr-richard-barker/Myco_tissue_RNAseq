# Original pre-analysis plan (superseded)

> **This is the plan written *before* the analysis was run, preserved for provenance.**
> It describes what was intended, not what happened. Several of its assumptions turned
> out to be wrong once measured — the aligner, the reference, the cause of low gene
> assignment, and the magnitude of the rRNA problem all changed. See [`../NOTES.md`](../NOTES.md)
> for every deviation with its evidence, and [`../README.md`](../README.md) for what was
> actually done.

---

Tissue-specific Pleurotus ostreatus RNA-seq: reanalysis, gene models, metabolic models
Context
16 RNA-seq libraries of Pleurotus ostreatus mycelial tissues were sequenced by Plasmidsaurus (project GPNJ7M). The vendor's entire delivered analysis is void: they aligned the reads to Saccharomyces cerevisiae. Proof is in GPNJ7M-expression-matrix.tsv — gene IDs are SGD systematic names (YDL151C/BUD30, YER174C/GRX4, Q0045/COX1), 7,127 features / 6,600 protein-coding = exactly SGD R64. Uniquely mapped reads are 453–2,284 per sample (0.006–0.03% of input); the ~16% "mapped" is multi-mapping noise into rRNA. The expression matrix, PCA, correlation heatmap, biotype plot and DGE results are all unusable. Their GSEA also used MSigDB Hallmark, which is human-only.

Independent confirmation the reads are fine: the most abundant read BLASTs at 100% identity to Pleurotus 28S LSU rRNA. The data are good; only the vendor's reference was wrong.

Goal: reprocess from raw FASTQ against a correct reference using a NASA GeneLab-derived pipeline, produce trustworthy per-locus counts and cross-sample normalisation, handle the heavy rRNA load with published justification, then build a per-tissue transcriptional model and a per-tissue context-specific metabolic model.

What the diagnostics established
Property	Value	Source
Assay	3'-end tag counting (QuantSeq-like), polyA-selected — not full-length RNA-seq	Plasmidsaurus RNA FAQ
Reads	SE 94 bp, NovaSeq X Plus, 14 nt UMI appended to read name after _	FASTQ headers
Batch	All 16 on one flowcell/run/lane (LH01080:180:25G52WLT4:4) — no batch term needed	FASTQ headers
Depth	2.28M–12.46M raw (median 5.4M); 11–62% of the vendor's advertised 20M	GPNJ7M-mapping-stats-reads.csv
rRNA	≥12–46%, varies systematically by tissue	4 diagnostic 25-mers over 200k reads/sample
Complexity	unique 50-mers only 7.5–35% of reads	same
Pre-processing	vendor already adapter-trimmed (0.002% residual TruSeq); 10.4% of reads carry ≥10 A	500k-read scan of 1A
Two consequences that force deviations from NF_RCP-F, detailed below: RSEM is the wrong estimator for 3'-tag data, and the PC9 RefSeq annotation has no UTRs (74,944 exons vs 74,754 CDS — exons ≈ CDS), so counting on exons alone would discard most 3'-tag reads.

Design
From Myco_Seq_Data_Plan.csv (corrected version; all 16 filenames verified against disk):

Tissue	Files	Wells	Raw reads (M)	rRNA floor (mean)
Exuding mycelium	GPNJ7M_1–_4	1A–1D	8.06, 5.58, 2.50, 5.47	31.5%
Fuzzy mycelium	GPNJ7M_5–_8	1E–1H	12.45, 4.92, 6.01, 5.38	36.1%
Exudophore	GPNJ7M_9–_12	2A–2D	2.74, 3.18, 3.80, 2.47	40.5%
Nodule	GPNJ7M_13–_16	2E–2H	2.28, 6.90, 3.34, 12.46	29.4%
One factor, 4 levels, n=4. Exudophore is the at-risk group: lowest depth and highest rRNA, so ~1.5M usable reads in the worst sample. rRNA fraction correlates with tissue, which is exactly why the rRNA-removed track is required rather than cosmetic — otherwise composition bias is confounded with the factor of interest.

Phase 0 — Project scaffold and environment
Create Tissue_specific_myeclium/analysis/ (git repo) with refs/, qc/, bam/, counts/, results/, models/, envs/, scripts/, plus metadata/runsheet.csv generated from Myco_Seq_Data_Plan.csv in GeneLab runsheet format.

Environment via miniforge/micromamba. Most bioconda RNA-seq tools have no osx-arm64 build, so create the alignment env with CONDA_SUBDIR=osx-64 and run under Rosetta 2:

env-align (osx-64): star=2.7.11b samtools=1.21 subread=2.1.1 umi_tools cutadapt=4.2 fastqc=0.12.1 multiqc rseqc=5.0.4 qualimap bedtools seqkit ucsc-gtftogenepred ucsc-genepredtobed
env-r (native arm64): r-base=4.4 bioconductor-deseq2=1.46 bioconductor-edger bioconductor-tximport r-tidyverse
env-model (native arm64): python=3.11 cobrapy modelseedpy troppo riptide memote
Version pins match GL-DPPD-7101-G so the work stays citable against the GeneLab spec. If Rosetta is unavailable, fall back to -profile docker for the alignment steps only.

Also: file a support ticket with Plasmidsaurus. Two grounds — wrong reference organism, and all 16 libraries delivered at 11–62% of the advertised 20M raw reads. Ask for a re-run against the correct genome and/or top-up sequencing, especially for the Exudophore group.

Phase 1 — Choose the reference empirically
No P. ostreatus GEM or obvious strain match exists a priori, and the LSU read confirms genus but not species. Decide with data. Subsample 500k reads/sample (seqkit sample), build minimal STAR indices, and map against each candidate:

Accession	Strain	Level	Size	Genes
GCF_014466165.1	PC9	contig (16)	34.9 Mb	11,849
GCA_029852705.2	PC9.15	chromosome (13)	35.3 Mb	13,556
GCA_056149245.1	BOM_ss5	contig (16)	41.4 Mb	12,705
GCA_056149315.1	BOM_ss14	contig (16)	43.5 Mb	13,310
GCF_014466165.1 is confirmed the RefSeq reference genome for the species (Academia Sinica, PC9, 2020) — your guess was right. But BOM_ss5/BOM_ss14 are Blue Oyster Mushroom single-spore isolates released March 2026 by the same group, and if your culture is a commercial blue oyster dikaryon they may match far better.

Select on: STAR % uniquely mapped + % reads assigned to genes (after Phase 2 annotation)

mismatch rate. Map to a single monokaryon, not a concatenation — for a dikaryon the two nuclear haplotypes are ~99% identical and concatenating them destroys unique mapping. If BOM_ss5 and BOM_ss14 score similarly and both beat PC9, that itself indicates a dikaryon; record it and pick the better single reference, noting allele-specific expression as out of scope. Also extract ITS from the reads to firm up the species call if the mapping rates are ambiguous.
Decision gate: if the best candidate still gives <60% unique mapping, stop and diagnose before proceeding — that would point to a strain not represented in public assemblies.

Phase 2 — Build a 3'-aware annotation and an rRNA reference
This is the technical crux and the main departure from the stock pipeline.

2a. rRNA reference. Run barrnap --kingdom euk on the chosen genome to locate 18S/5.8S/28S/5S, add the mitochondrial rRNAs (MT is present, 73 kb), and supplement with an rDNA unit assembled from the reads themselves (rDNA repeats are usually collapsed or absent in assemblies, which is why the vendor's rRNA reads had nowhere legitimate to go). Emit refs/rRNA_loci.bed plus a FASTA for direct read-level rRNA quantification. This gives the honest rRNA fraction, replacing the 4-probe floor.

2b. PAS atlas and 3'-end extension. Because the annotation has zero UTRs:

Align once with soft-clipping retained (Phase 3 settings).
Call polyadenylation sites from reads whose soft-clipped tail is untemplated A — these pinpoint the actual cleavage site. 10.4% of reads carry ≥10 A, so there is ample signal.
Filter internal priming: drop any PAS whose downstream genomic 20 bp window has ≥6 consecutive A or ≥12/20 A. Standard and essential for polyA-primed 3' methods.
Assign surviving PAS to the nearest upstream gene on the same strand.
Extend each gene's 3' end to its furthest supported PAS, capped at min(500 bp, distance to next gene on the same strand − 50 bp). Fungal genes are compact, so the neighbour cap matters — without it you get read-through cross-assignment.
Genes with no PAS support get a flat +300 bp extension under the same cap.
Emit refs/annotation_3p.gtf. Acceptance check: fraction of reads assigned by featureCounts before vs after extension. Expect a large jump; if it does not rise materially, the extension logic is wrong and must be fixed before any counts are trusted.

Side benefit: this PAS atlas is, as far as the literature search shows, the first for P. ostreatus, and it enables alternative-polyadenylation comparison between tissues — a genuinely novel result available at no extra cost.

Phase 3 — Adapted GeneLab pipeline (per GL-DPPD-7101-G)
Note NF_RCP-F is superseded by GL-DPPD-7101-G (Feb 2025), which added exactly the parallel rRNA-removed DGE track you asked about. Follow G, not F.

Raw QC — FastQC 0.12.1 + MultiQC, confirm the vendor's trimming state.
Trim — cutadapt: polyA trimming (--poly-a), quality trim, -m 30. Light adapter pass only; do not re-trim aggressively over already-trimmed reads. Keep the untrimmed BAM for PAS calling (Phase 2b needs the soft-clipped A tails).
STAR index — --genomeSAindexNbases 11 (35 Mb genome), --sjdbGTFfile, --sjdbOverhang 93.
STAR align (SE) — --outFilterMultimapNmax 20 --alignIntronMax 3000 (fungal introns are short; the default lets reads span implausible gaps), --outSAMattributes NH HI AS nM, --outReadsUnmapped Fastx (the G-version addition).
UMI dedup — umi_tools dedup --umi-separator=_ --method=directional. The 14 nt UMI is in the read name. Position collisions are the norm for 3'-tag data, so directional error-aware collapsing matters.
Post-align QC — RSeQC infer_experiment (settles strandedness empirically), read_distribution, geneBody_coverage (expect a strong 3' skew — that is correct here, not a failure), Qualimap; BED built via gtfToGenePred/genePredToBed as GeneLab does.
Quantify — featureCounts -t exon -g gene_id -s <inferred> against refs/annotation_3p.gtf. Primary counts use uniquely-mapping reads; a -M --fraction sensitivity run is reported alongside.
Phase 4 — Normalisation and differential expression
No TPM/FPKM. 3'-tag counting yields one count per molecule independent of transcript length, so effective-length normalisation is meaningless. Counts and CPM only. This alone invalidates any length-normalised summary of these data.
Two parallel tracks, mirroring GL-DPPD-7101-G: all-genes and rRNA-removed (*_rRNArm), each independently re-normalised after filtering. Justification to state in the manuscript: rRNA is 12–46% of reads and varies with tissue (Exudophore ~40% vs Nodule ~29%), so retaining it imposes a composition bias confounded with the design factor; NASA GeneLab's current consensus pipeline performs this same removal by default.
edgeR::filterByExpr, then DESeq2 median-of-ratios; VST counts for PCA/clustering (the G-version addition).
Model ~ Tissue; all 6 pairwise contrasts plus an LRT for any-tissue effect.
Diagnostics: library size vs PC1, rRNA fraction vs PC1, dedup rate per sample. Given the 2.3–12.5M depth spread, confirm no PC tracks depth rather than biology.
Sanity gate before interpretation: replicates must cluster by tissue on VST-PCA. If Exudothore replicates scatter, its low usable depth is the likely cause and that group's conclusions get an explicit power caveat.

Phase 5 — Per-tissue transcriptional model
Interpreting "genetic model" as a tissue-resolved expression/regulatory model — flag if you meant something narrower.

Functional annotation of the chosen proteome: eggNOG-mapper (KO, GO, COG), InterProScan, dbCAN3 for CAZymes (central to a wood-decay fungus), SignalP + DeepTMHMM for the secretome, and manual curation of laccase / Mn-peroxidase / versatile-peroxidase families.
Per-tissue expressed gene set at a CPM threshold; tau specificity index; tissue-enriched and tissue-specific classification.
Marker gene sets per tissue with effect sizes.
Co-expression modules on VST counts. Caveat: n=16 is marginal for WGCNA — report modules as exploratory, or use a method better behaved at low n.
Enrichment per tissue against the eggNOG-derived KEGG/GO sets (not Hallmark).
Expected biology to test rather than assume: exuding mycelium and exudophore should show elevated secretome, transporters, and polyol/osmolyte metabolism (mannitol, trehalose, oxalate) relative to fuzzy mycelium and nodule.

Phase 6 — Per-tissue metabolic models
No published P. ostreatus GEM was found, so this is a de novo reconstruction — a publishable contribution in its own right.

Draft reconstruction — ModelSEEDpy with a fungal template from the annotated proteome (pip-installable, runs natively on the Mac; KBase's web build is the no-install fallback).
Compartmentalisation — DeepLoc2 / WoLF PSORT for cytosol / mitochondrion / peroxisome assignment. Non-optional for fungal metabolism, where β-oxidation and glyoxylate cycle placement change flux predictions substantially.
Gapfill on a defined medium matching the actual culture substrate — needed from you.
QC with memote; curate the lignocellulose and secondary-metabolism blocks by hand, since automated templates cover them poorly.
Context-specific extraction per tissue with RIPTiDe (transcriptome-weighted flux sampling, well suited to 4 conditions) cross-checked against Troppo (GIMME/iMAT/tINIT). Input: rRNA-removed normalised counts.
Analysis — FBA/FVA per tissue; reaction presence/absence contrasts; and predicted secretion capacity, which ties the exudophore and exuding-mycelium models directly to droplet composition.
Deviations from NF_RCP-F, and why
GeneLab default	Here	Reason
NF_RCP-F / GL-DPPD-7101-F	GL-DPPD-7101-G	G is current and adds the rRNA-removed DGE track
Ensembl reference auto-fetch	NCBI RefSeq/GenBank, chosen empirically	P. ostreatus is not in Ensembl
RSEM quantification	featureCounts on a 3'-extended annotation	RSEM models full-length transcript coverage; invalid for 3'-tag pileups
Annotation used as-is	PAS-derived 3' extension	RefSeq PC9 GTF has no UTRs; 3'-tag reads land beyond the stop codon
TPM/FPKM reported	counts/CPM only	no effective-length normalisation for 3'-tag data
Paired-end handling	single-end throughout	SE 94 bp library
Nextflow + Singularity	conda envs, version-pinned to G	Singularity has no native macOS support
Trim Galore!	cutadapt with polyA trimming	reads are pre-trimmed; polyA is the actual problem
(no UMI step)	umi_tools directional dedup	14 nt UMI present and essential for 3'-tag
Verification
Reference choice — mapping-rate table across the 4 candidates; the winner must exceed 60% unique.
Annotation fix — featureCounts assignment rate before vs after 3' extension, and the PAS internal-priming filter pass rate.
rRNA — read-level rRNA fraction from the Phase 2a reference, replacing the probe-based floor.
Counts — RSeQC read_distribution showing exonic/3'-biased placement; per-sample assignment rate; total counts vs input reads reconciled.
DE — VST-PCA clusters by tissue; replicate correlation within group exceeds between; all-genes and rRNA-removed tracks agree in direction for top genes.
Models — memote score on the draft GEM; each tissue model must produce non-zero biomass flux on the defined medium; marker genes from Phase 5 must appear in the corresponding tissue model's active reaction set.
End-to-end — one MultiQC report covering raw, trimmed, aligned and dedup stages.
Deliverables
Reproducible analysis/ repo; corrected expression matrix (counts, CPM, VST) for all 11.8k–13.5k loci; rRNA-removed parallel matrix; DE tables for 6 contrasts; per-tissue marker and specificity tables; functional annotation of the proteome; a P. ostreatus PAS atlas; a draft GEM plus 4 tissue-specific models; and a methods section written against GL-DPPD-7101-G with every deviation justified.

Open items — needed from you
Culture substrate / medium composition — required to gapfill and constrain the GEMs (Phase 6.3).
Strain provenance — commercial blue oyster isolate, or a lab strain? Informs Phase 1 even though the decision is made on mapping rates.
Confirm "genetic model" means the tissue-resolved expression/regulatory model of Phase 5.
Risks
Underpowered Exudophore group. Lowest depth, highest rRNA, ~1.5M usable reads at worst. Top-up sequencing from the vendor is the clean fix; otherwise this group carries a power caveat.
No public strain match. If Phase 1 fails its 60% gate, a de novo transcriptome assembly becomes the fallback — a materially larger scope change, which I would bring back to you.
n=4 with one factor limits co-expression network inference; treated as exploratory.
