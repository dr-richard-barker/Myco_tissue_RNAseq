# Ribosomal RNA dominance and unannotated mitochondrial rDNA confound 3′-tag RNA-seq in *Pleurotus ostreatus*: a corrected reference, annotation and metabolic reconstruction

**Draft — resource/methods paper. Target: fungal-biology venue** (e.g. *Fungal Genetics and
Biology*, *Fungal Biology and Biotechnology*, or *IMA Fungus*).

> **Status of this draft.** All numbers are generated from the analysis in this repository and
> are traceable to the files named in each section. Nothing here is estimated or carried over
> from planning documents. Remaining **[TO WRITE]** markers indicate the few items only the
> authors can supply — chiefly the sampling protocol and the morphological description of the
> exudophore. The tissue
> biology is deliberately limited to a demonstration of pipeline function; the full
> tissue-comparison analysis is reserved for a companion paper.

---

## Abstract

Transcriptomic studies of non-model basidiomycetes increasingly rely on 3'-end tag RNA-seq and
on public genome assemblies whose annotations were built for gene finding rather than for read
quantification. We show that both choices carry substantial and largely silent costs. From
92.4 million reads across 16 libraries of *Pleurotus ostreatus* cv. "Harbor Blue P01" grown in
a mycoponic ceramic-tube system, only 2.48 million reads (2.7%) were assignable to
protein-coding features. Nuclear ribosomal DNA accounted for a mean of 53.5% of reads, and the
single largest remaining component, mitochondrial rRNA, **is unannotated in all four public
*P. ostreatus* genome annotations we examined**. Because it carries no feature it is neither
counted nor removable by the biotype-based rRNA filtering that current consensus pipelines
implement; it is silently discarded as "no feature". We further show that the intuitive remedy
for low assignment in tag data — extending gene 3' ends — recovered 0.1% here, and we give the
diagnostic that distinguishes truncated UTRs from unannotated rRNA. Assembly choice also
proved consequential in an unexpected way: the RefSeq-designated reference genome for the
species performed worst of four candidates (5.45% versus 24.94% uniquely mapped), and
assemblies differ in whether they resolve or collapse the nuclear rDNA array, which changes
multi-mapping and therefore quantification. We provide a corrected, rRNA-complete and 3'-aware
annotation, a functional annotation of the proteome including 423 CAZymes and a predicted
secretome, and the first draft genome-scale metabolic reconstruction for *P. ostreatus*, which
achieves non-zero biomass flux on the defined culture medium. Despite the compromised input,
the corrected pipeline recovered convergent biological signal in two tissues by two
independent criteria, demonstrating both what such data can still support and where the limit
lies.

---

## 1. Introduction

*Pleurotus ostreatus* is among the most widely cultivated edible fungi and a workhorse model
for white-rot lignocellulose degradation. Its enzymatic repertoire — laccases, manganese and
versatile peroxidases, lytic polysaccharide monooxygenases and an extensive glycoside
hydrolase complement — underpins both its ecological role and a growing set of biotechnological
applications. Interest has recently extended to controlled-environment and bioregenerative life
support contexts, where fungal biomass offers a protein source that can be produced from
lignocellulosic residues in a closed system. The mycoponic culture format used here, in which
mycelium is grown on micro-structured ceramic tubes supplied with liquid nutrient medium
through an antimicrobial size-exclusion interface, was developed for exactly this purpose
(Porterfield et al. 2026).

Understanding how such a system works requires resolving what different parts of the mycelium
are doing. Filamentous fungal colonies are not homogeneous: they differentiate into
morphologically and physiologically distinct regions, and in this culture system they form
several visually distinguishable tissue types, including structures that produce and exude
liquid droplets. Transcriptome comparison across such tissues is the natural first approach.

For multi-condition designs of this kind, 3'-end tag counting has become an attractive option:
it sequences a single fragment per transcript, so cost per sample falls sharply and library
complexity requirements are lower than for full-length protocols. The trade-offs are less
widely appreciated than the benefits. Because one read corresponds to one molecule regardless
of transcript length, no length normalisation applies and TPM or FPKM values are meaningless.
Because reads pile up at the polyadenylation site, quantification depends on the accuracy of
annotated 3' ends rather than on gene bodies. And because the method primes on poly(A) tracts,
it is sensitive to internal priming and to any failure of poly(A) selection.

Those dependencies interact badly with a second issue: the state of public fungal genome
annotation. Genome annotations are produced to describe protein-coding gene models. Ribosomal
DNA is frequently omitted, being repetitive, hard to assemble and of little interest for gene
finding; organellar rRNA is omitted more often still. For a genome browser this is
inconsequential. For read quantification it is not, because a feature that does not exist
cannot receive counts and cannot be filtered out. Current consensus RNA-seq pipelines,
including the NASA GeneLab pipeline used as the basis for this work, implement rRNA removal by
dropping features whose annotated biotype is rRNA — a strategy that is silently ineffective
against rRNA that was never annotated.

We encountered both problems in an acute form. The libraries analysed here were commercially
generated, delivered with an analysis aligned to the wrong species, and proved on reanalysis
to be dominated by ribosomal RNA. Rather than treat this as a private failure, we use it to
quantify costs that are usually invisible, and to provide the corrected resources that the
species has lacked. Specifically, we ask: how much does reference choice matter for a species
with several published assemblies; how much signal is lost to unannotated rRNA and how can it
be detected; is the intuitive explanation for poor assignment in tag data the correct one; and
what can and cannot be concluded from data compromised in this way.

## 2. Materials and methods

### 2.1 Biological material and culture
*Pleurotus ostreatus* cv. "Harbor Blue P01" was grown on micro-structured ceramic MiniTubes
(10 cm × 5 cm, pore size < 300 nm, 50% v/v granular activated carbon) supplied with mycoponic
nutrient medium v3, at 16 °C, 85% relative humidity and CO₂ ≤ 1000 ppm, as described by
Porterfield et al. (2026). Four tissue types were sampled with four biological replicates
each: exuding mycelium, fuzzy mycelium, **exudophore** (a newly described exudate-producing
structure), and nodule. Full medium composition and culture parameters are recorded in
`metadata/culture_conditions.md`.

**[TO WRITE — authors]** Precise sampling procedure, tissue definitions and morphological
description of the exudophore, RNA extraction method, and age of cultures at sampling.

### 2.2 Sequencing, and its quality

Libraries were prepared and sequenced commercially as poly(A)-selected 3'-end tag libraries
(Illumina NovaSeq X Plus, single-end 94 bp, 14 nt unique molecular identifier appended to the
read name). All 16 libraries were sequenced on one flowcell and lane
(`LH01080:180:25G52WLT4:4`), so no batch term was required. Total yield was 92,430,743 raw
reads (median 5.4 M per library).

We report the following quality observations in full, because they materially constrain the
conclusions that follow and because they are not evident from the delivered quality-control
outputs.

**The delivered analysis was aligned to the wrong species.** The provider's expression matrix
contains *Saccharomyces cerevisiae* systematic gene identifiers (`YDL151C`, `YER174C`,
`Q0045`/COX1); its 7,127 features and 6,600 protein-coding genes correspond exactly to the SGD
R64 annotation. Uniquely mapped reads numbered 453–2,284 per library, or 0.006–0.03% of input.
The accompanying gene-set enrichment analysis used MSigDB Hallmark, a human collection. The
delivered expression matrix, principal component analysis, correlation heatmap, biotype
summary and differential expression results are consequently without meaning. The reads
themselves are sound: the most abundant read matches *Pleurotus* 28S rRNA at 100% identity.

**Sequencing depth fell short of specification.** All 16 libraries were delivered at 11–62% of
the 20 M raw reads per sample advertised for the service (median 27%).

**Poly(A) selection did not enrich mRNA.** Ribosomal RNA constituted 20–70% of reads per
library (Section 3.2). This is a property of the delivered libraries rather than of the
downstream analysis, and it is the dominant constraint on the dataset.

These observations are separable in kind. The species mis-assignment is an analytical error in
the delivered report and does not affect the underlying reads. The depth shortfall is a
quantitative deviation from specification. The ribosomal content is a library-preparation
outcome that neither party detected before delivery, and which the provider's own summary
statistics — reporting approximately 16% of reads as "mapped" without comment — would have
flagged had the reference been correct. We note these distinctions because they carry different
implications, and because the general lesson for users of commercial sequencing is that
provider quality-control metrics computed against a wrong reference are not merely uninformative
but actively misleading.

### 2.3 Reference selection
Four *P. ostreatus* assemblies were compared empirically rather than assumed
(`scripts/02_testmap.sh`): PC9 (`GCF_014466165.1`), PC9.15 (`GCA_029852705.2`), BOM_ss5
(`GCA_056149245.1`) and BOM_ss14 (`GCA_056149315.1`). 250,000 reads per library were mapped to
genome-only indices and compared on uniquely mapped fraction, then on gene assignment rate
using matched 3′-extended annotations.

### 2.4 Annotation construction
Ribosomal RNA loci were identified with barrnap (nuclear and mitochondrial modes) and
supplemented with coverage-derived blocks obtained by clustering intergenic reads from pooled
alignments (`scripts/04_rrna_ref.sh`, `refs/rRNA_regions_<label>.tsv`). Gene 3′ ends were
extended strand-aware by up to 500 bp, capped at the distance to the nearest neighbouring
feature minus 50 bp (`scripts/07_extend_3p.py`). rRNA blocks were added as `gene_biotype
"rRNA"` features on a single strand (`scripts/08_build_annotation.py`).

### 2.5 Read processing and quantification
Adapted from the NASA GeneLab RNAseq Consensus Pipeline GL-DPPD-7101-G. Reads were trimmed
with fastp 1.3.6 (polyG, polyX, 3′ quality, minimum length 30), aligned with HISAT2 2.2.3
(`--max-intronlen 3000 --rna-strandness F`), reduced to primary alignments, deduplicated with
umi_tools 1.1.6 (`--method=directional`), and counted with featureCounts (subread 2.1.1;
`-t exon -g gene_id -s 1`). Strandedness was established empirically (forward: 49.6% assigned
versus 6.3% reverse).

**Deviations from GL-DPPD-7101-G**, each with its rationale, are tabulated in `NOTES.md` §7.
The two consequential ones: HISAT2 replaces STAR because STAR is non-functional on the
analysis platform (Darwin 25.5 / Apple Silicon; it loads a genome then reports zero input
reads, reproducibly, including on a synthetic control); and featureCounts replaces RSEM
because RSEM models full-length transcript coverage and is invalid for 3′-tag pileups. No
TPM/FPKM values are reported anywhere: 3′-tag counting yields one count per molecule
independent of transcript length.

### 2.6 Differential expression
DESeq2 1.46.0 with `~ Tissue`, following GL-DPPD-7101-G in running two parallel tracks
(all-genes and rRNA-removed, each independently re-normalised). Genes were filtered with
`edgeR::filterByExpr`. Libraries yielding fewer than 15,000 assigned non-rRNA counts were
excluded in a sensitivity analysis (Section 3.5).

### 2.7 Functional annotation
DIAMOND blastp against UniProt Swiss-Prot; EC, GO (with ontology aspect) and KEGG pathway
assignments derived from the matched Swiss-Prot records, with EC→pathway mapping via the KEGG
REST API. CAZy classes were assigned from Pfam-A domains using curated keyword rules. Secreted
proteins were inferred by transferring curated SIGNAL and TRANSMEM features from the best-hit
Swiss-Prot entry (signal peptide present, no transmembrane helix beyond it).

**Substitutions and their limitations.** dbCAN was unavailable (its download service redirects
to a landing page; no database file was retrievable), so CAZy assignments are **Pfam-derived
and class-level only**, not dbCAN family calls. SignalP and Phobius are licence-restricted and
DeepSig's distribution pins an unavailable TensorFlow version, so the secretome is
**homology-transferred rather than predicted de novo**, covering only proteins with a
Swiss-Prot hit and therefore representing a lower bound. eggNOG-mapper's database host was
unreachable, so Swiss-Prot was used instead.

### 2.8 Co-expression network
WGCNA 1.74, signed network, soft-threshold power 18, minimum module size 30, merge height
0.25. Module–trait correlations were corrected across all module × tissue tests using
Benjamini–Hochberg. Enrichment of GO (BP/MF/CC) and KEGG terms per module used the
hypergeometric test with BH correction within each module × ontology.

### 2.9 Metabolic reconstruction
EC numbers from the functional annotation were mapped onto the ModelSEED biochemistry database
to produce a draft reconstruction with gene–protein–reaction rules traceable to supporting
proteins. Transport reactions were added and exchange bounds constrained to the mycoponic
nutrient medium. Analyses used COBRApy 0.32.1.

### 2.10 Availability
All code, count matrices, annotation, models and intermediate results are in the repository;
see `README.md`, `metadata/data_dictionary.md` and `metadata/reference_sources.tsv`. Conda
environments are pinned in `envs/locks/`. **[TO ADD]** Raw reads: SRA/ENA accession pending.
Archived release: Zenodo DOI pending.

---

## 3. Results

### 3.1 Public annotations for this species mis-assign the reference by an order of magnitude

Uniquely mapped fraction across the four candidate assemblies (250k reads × 16 libraries,
genome-only indices) differed far more than assembly statistics would suggest:

| Assembly | Strain | Level | Unique % | Multi % | Unaligned % |
|---|---|---|---|---|---|
| BOM_ss5 (`GCA_056149245.1`) | blue oyster SSI | contig | **24.94** | 53.17 | 21.89 |
| BOM_ss14 (`GCA_056149315.1`) | blue oyster SSI | contig | 24.92 | 53.18 | 21.90 |
| PC9.15 (`GCA_029852705.2`) | PC9.15 | chromosome | 21.49 | 53.10 | 25.42 |
| PC9 (`GCF_014466165.1`) | PC9 | contig | **5.45** | 69.14 | 25.41 |

The RefSeq-designated reference genome for the species (PC9) performed **worst**, at less than
a quarter the unique mapping rate of the best candidate. BOM_ss5 and BOM_ss14 differed by 0.02
percentage points, consistent with their being the two nuclei of a single dikaryon — matching
the cultivar used here.

On matched 3′-extended annotations, BOM_ss5 assigned 6.9% more reads than PC9.15 and did so
for **every one of the 16 libraries** (range +4.0% to +13.5%). BOM_ss5 was used for all
subsequent analysis; the complete parallel PC9.15 analysis is retained in the repository.

### 3.2 The libraries are dominated by rRNA from two distinct sources

Against barrnap-derived loci, nuclear rDNA accounted for a mean of 53.5% of reads
(range 19.9–69.7%). This is a direct measurement; a k-mer probe estimate made earlier in the
analysis (12–46%) was a substantial underestimate, and we note this because probe-based
estimates are common in QC reporting.

Clustering intergenic reads from pooled alignments revealed the second source. In BOM_ss5,
**40.2% of pooled aligned reads fall in a tandem block structure at
`JBQVBD010000012.1:2,298,112–2,346,320`** (the nuclear rDNA array) and **31.3% on the
mitochondrial contig `CM148777.1`** — approximately 71% of aligned reads in identifiable rDNA
loci. barrnap's fragmentary 5S/5.8S/18S/28S calls fall inside every block, confirming their
identity.

In PC9.15 the equivalent mitochondrial block (`CM057219.1:52,604–56,069`) contained **60.9% of
all pooled aligned reads** and was confirmed by BLAST at 100% identity to the *P. ostreatus*
mitochondrion.

**Neither the nuclear rDNA array nor the mitochondrial rRNA is annotated as a feature in any of
the four public annotations.** In the GTF, the mitochondrial rRNA region is flanked only by
tRNA gene models. The practical consequence is severe: these reads cannot be counted, cannot be
removed by the biotype-based rRNA filtering that current consensus pipelines implement, and
are silently reported as "no feature".

### 3.3 Assembly choice changes rDNA behaviour, not just contiguity

BOM_ss5 resolves multiple tandem copies of the nuclear rDNA repeat; PC9.15 collapses them. The
consequence is quantitative, not cosmetic: in the collapsed assembly, rDNA reads distribute
across fewer loci and are reported as multi-mapping (and discarded), whereas in the resolved
assembly they pile into an identifiable block that can be annotated and accounted for. This is
the mechanism underlying the assignment-rate difference in Section 3.1, and it means that
**assembly selection for RNA-seq should consider rDNA representation, which is not reported in
standard assembly metrics.**

### 3.4 Extending gene 3′ ends does not fix the problem, and the diagnostic that reveals why

A natural response to low assignment in 3′-tag data is that annotated 3′ UTRs are too short.
This is superficially well supported here: of the four annotations, only PC9.15 carries
meaningful UTRs (median 3′ UTR 78 bp; the other three have 3 bp, i.e. the stop codon only).

We extended gene 3′ ends strand-aware by up to 500 bp with a neighbour cap (91.3% of 13,556
genes extended; median extension 323 bp; 3.79 Mb added). **The effect on gene assignment was
15.5% → 15.6%.**

The diagnostic that distinguishes the two causes is simple and we recommend it as routine:
**cluster the reads that fail to be assigned and ask whether they are diffuse or focal.** UTR
truncation produces diffuse loss distributed just downstream of many genes; unannotated rRNA
produces a small number of very high-coverage blocks. Here a single locus held 60.9% of
aligned reads. Adding rRNA features reduced "no feature" from 30–78% to 0–1.2%.

### 3.5 Effective yield, and the cost of retaining noise-dominated libraries

Despite 92.4 M raw reads, only 2,475,808 reads were assigned to non-rRNA features across all 16
libraries; 7,601,974 were assigned to rRNA features. Effective mRNA yield ranged from
**0.09% to 10.49% of raw reads**, and only 3 of 16 libraries exceeded 100,000 assigned
non-rRNA counts.

**Table 1.** Per-library read budget (BOM_ss5 reference, UMI-deduplicated).

| Tissue | Well | Raw reads | rRNA counts | mRNA counts | mRNA % raw | Genes >=10 | Retained |
|---|---|---|---|---|---|---|---|
| Exuding mycelium | 1A | 8,423,109 | 701,579 | 469,981 | 5.58 | 6,017 | yes |
| Exuding mycelium | 1B | 5,902,765 | 459,570 | 19,365 | 0.33 | 354 | yes |
| Exuding mycelium | 1C | 2,770,639 | 64,744 | 2,620 | 0.09 | 19 | **no** |
| Exuding mycelium | 1D | 5,770,952 | 646,755 | 29,276 | 0.51 | 546 | yes |
| Exudophore | 2A | 2,941,102 | 145,479 | 13,447 | 0.46 | 239 | **no** |
| Exudophore | 2B | 3,412,844 | 99,368 | 10,944 | 0.32 | 164 | **no** |
| Exudophore | 2C | 4,028,384 | 143,150 | 66,739 | 1.66 | 1,485 | yes |
| Exudophore | 2D | 2,620,662 | 168,385 | 46,606 | 1.78 | 998 | yes |
| Fuzzy mycelium | 1E | 13,369,308 | 1,032,480 | 35,940 | 0.27 | 704 | yes |
| Fuzzy mycelium | 1F | 5,289,692 | 330,615 | 46,085 | 0.87 | 1,034 | yes |
| Fuzzy mycelium | 1G | 6,274,968 | 714,848 | 224,213 | 3.57 | 4,507 | yes |
| Fuzzy mycelium | 1H | 5,610,507 | 887,407 | 81,814 | 1.46 | 1,698 | yes |
| Nodule | 2E | 2,431,986 | 118,075 | 7,062 | 0.29 | 109 | **no** |
| Nodule | 2F | 7,338,750 | 588,323 | 41,959 | 0.57 | 827 | yes |
| Nodule | 2G | 3,473,216 | 531,658 | 40,253 | 1.16 | 813 | yes |
| Nodule | 2H | 12,771,859 | 969,538 | 1,339,504 | 10.49 | 7,512 | yes |


Retaining the lowest-yield libraries actively degraded the analysis. With all 16, the first
principal component correlated with sequencing depth at r = 0.686 — that is, the dominant axis
of variation was library size, not biology. Excluding the four libraries below 15,000 assigned
non-rRNA counts (1C, 2E, 2B, 2A) reduced this to **r = −0.003**, while increasing genes passing
expression filtering from 1,079 to 1,666 and genes with a significant tissue effect from 76 to
281. All four tissues retained n ≥ 2.

### 3.6 Resources generated

**Annotation.** rRNA-complete, 3′-aware GTFs for BOM_ss5 and PC9.15
(`refs/<label>/<label>_final.gtf`), with rRNA blocks declared in per-reference TSVs.

**Functional annotation.** Of 12,521 BOM_ss5 proteins, 5,163 (41%) have a Swiss-Prot hit;
39,607 GO/KEGG term assignments; **423 proteins carrying a CAZy-class domain** (GH 233,
**AA 86**, GT 57, PL 30, CBM 13, CE 4); **576 proteins predicted secreted**. The AA (auxiliary
activity) complement — laccases, peroxidases, LPMOs — is directly relevant to lignocellulose
degradation and is consistent with the oak sawdust and microcrystalline cellulose in the
growth medium.

**Metabolic reconstruction.** The first draft genome-scale metabolic model for
*P. ostreatus*: 1,915 reactions, 2,028 metabolites and 881 genes, with GPR rules traceable to
supporting proteins (`models/BOM_ss5_draft.xml`). Constraining exchanges to the mycoponic
medium and adding transport reactions yields 5,242 reactions of which **2,262 can carry flux**
(from 288 in the unconstrained draft), and 28 of 34 standard biomass precursors become
synthesisable, including all 20 proteinogenic amino acids and all NTPs/dNTPs
(`models/BOM_ss5_medium.xml`).

Targeted gapfilling then completed the cofactor pathways. Restricting the candidate pool by
backward reachability from each target metabolite — rather than searching the whole ModelSEED
universal database — made each target solvable in seconds, and required only **three added
reactions with no gene evidence** (tagged as gapfilled). NAD, FAD and UDP-N-acetylglucosamine
each needed one reaction; CoA and NADP were already producible. Biotin remained unreachable
and was added to the medium rather than given a biosynthesis route, which is biologically
appropriate: many fungi are biotin auxotrophs and the medium's malt extract, peptone and
tryptic soy broth all supply B-vitamins. The resulting model
(`models/BOM_ss5_gapfilled.xml`) has 5,247 reactions of which 2,287 carry flux, all six tested
cofactors producible, and **non-zero biomass flux on the defined medium**.

**This model remains a scaffold, not a validated predictor.** The biomass objective is coarse
and uncurated, mass balance has not been curated, and three reactions are present on
topological rather than genomic evidence.

### 3.7 Demonstration: the pipeline recovers coherent signal where power permits

To show the corrected pipeline yields interpretable biology, we assessed tissue signal by two
independent criteria.

**Replicate-supported markers.** Tissue-specificity was scored with the τ index on linear
normalised counts. Because τ is computed on tissue means, a single outlying replicate can
manufacture an apparent marker; we therefore required that **every** replicate of a tissue
exceed the maximum value in any sample of any other tissue. Of the top 50 markers per tissue:

| Tissue | n | Markers supported by all replicates |
|---|---|---|
| Exudophore | 2 | **37 / 50** |
| Nodule | 3 | 22 / 50 |
| Exuding mycelium | 3 | 1 / 50 |
| Fuzzy mycelium | 4 | 1 / 50 |

Note that the test becomes more stringent as n increases, so these percentages are not
comparable across tissues; the qualitative separation is nonetheless clear.

**Co-expression modules.** WGCNA (12 libraries, 4,748 genes) produced 41 modules, 37 with at
least one enriched GO or KEGG term. Of 168 module × tissue correlations, 17 reached p < 0.05
(8.4 expected by chance) and **6 survived BH correction — all associated with exudophore or
nodule**:

| Module | Genes | Enriched terms (title) | Tissue | r | FDR |
|---|---|---|---|---|---|
| brown | 320 | oxidoreductase activity, CH-OH donors; extracellular space | Nodule | +0.85 | 0.016 |
| lightcyan | 99 | SSU-rRNA maturation; CENP-A chromatin assembly | Exudophore | +0.91 | 0.003 |
| lightsteelblue1 | 34 | chromosome segregation | Exudophore | +0.85 | 0.016 |
| saddlebrown | 49 | methane metabolism; microbial metabolism in diverse environments | Exudophore | +0.84 | 0.018 |
| greenyellow | 142 | pexophagy; response to calcium ion | Exudophore | −0.92 | 0.003 |
| white | 55 | one-carbon metabolic process | Nodule | −0.88 | 0.009 |

The two methods agree completely on which tissues carry signal, despite being sensitive to
different failure modes. Neither supports exuding or fuzzy mycelium at the achieved depth.

**Convergent signatures.** The nodule-associated `brown` module enriches for extracellular
space (7/20 genes, FDR = 0.003) and CH-OH-acting oxidoreductases, and contains 22 predicted
secreted proteins and 11 CAZymes; it also enriches for "protein aggregate center" (3/5). The
replicate-supported nodule markers are independently dominated by cerato-platanin family
proteins and small heat-shock proteins — consistent with the aggregate-centre enrichment
recovered by a different method.

The **exudophore**, a newly described exudate-producing structure, yields the largest set of
replicate-supported markers of any tissue. These are consistently oxidative and secretory:
alcohol oxidase 1, glyoxal oxidase, formate and aldehyde dehydrogenases, an FAD-dependent
monooxygenase, acetate–CoA ligase, an amino-acid permease and an AA9 LPMO. Alcohol oxidase and
glyoxal oxidase both generate extracellular H₂O₂, the co-substrate required by fungal
peroxidases.

### 3.8 Tissue modules are not explained by ribosomal RNA load

Because rRNA content varies systematically among these libraries (42.0–96.6% of assigned
counts) and is itself correlated with usable depth (r = −0.907 with log₁₀ mRNA counts), any
module correlating with tissue could in principle reflect rRNA load rather than biology. This
concern is sharpest for the exudophore-associated `lightcyan` module, which enriches for
SSU-rRNA maturation — a ribosome-biogenesis signature that could plausibly co-vary with
ribosomal content for purely technical reasons.

We tested each module eigengene against the per-sample rRNA fraction and recomputed the
tissue association as a partial correlation controlling for it
(`scripts/25_rrna_confound.R`):

| Module | Tissue | r (tissue) | r (rRNA fraction) | r (partial) | FDR (partial) |
|---|---|---|---|---|---|
| greenyellow | Exudophore | −0.92 | +0.28 | −0.92 | 0.006 |
| lightcyan | Exudophore | +0.91 | **−0.15** | +0.91 | 0.006 |
| white | Nodule | −0.88 | −0.10 | −0.92 | 0.006 |
| lightsteelblue1 | Exudophore | +0.85 | −0.36 | +0.85 | 0.025 |
| brown | Nodule | +0.85 | −0.09 | +0.85 | 0.025 |
| saddlebrown | Exudophore | +0.84 | +0.07 | +0.88 | 0.015 |

All six tissue associations survive. The `lightcyan` ribosome-biogenesis module correlates
with rRNA fraction at only −0.15, so there is essentially no technical signal to remove, and
its exudophore association is unchanged by the correction. Because rRNA fraction is strongly
anti-correlated with usable depth, this control also partially accounts for depth.

**We report these as evidence that the corrected pipeline recovers coherent, method-independent
signal, not as a characterisation of exudophore function.** With n = 2 after quality filtering,
the full tissue comparison is reserved for a companion study.

---

## 4. Discussion

### 4.1 Unannotated rRNA is a silent failure mode, not a noisy one

The central practical finding is that the largest single component of these libraries occupied
a genomic region carrying no annotated feature. Its consequences were entirely silent. No tool
raised an error. Alignment rates looked unremarkable. The reads were reported by featureCounts
under `Unassigned_NoFeatures`, a category most workflows summarise but few interrogate, and
which is easily rationalised as intergenic transcription or annotation incompleteness in the
diffuse sense.

The consequence for rRNA removal deserves emphasis. Version G of the NASA GeneLab consensus
pipeline added a parallel rRNA-removed differential expression track, implemented by dropping
features whose biotype is rRNA before renormalisation. That is a sound design, and we adopted
it. But it is a *feature*-based operation, and it therefore cannot remove rRNA that was never
annotated as a feature. A pipeline can execute its rRNA-removal step faithfully, report success,
and leave the majority of the ribosomal signal untouched. Any workflow implementing
biotype-based rRNA filtering inherits this vulnerability, and the vulnerability grows precisely
in the situation where it matters most — a non-model organism whose annotation was contributed
by a genome project rather than curated for expression analysis.

Two properties made this detectable. First, the loss was focal rather than diffuse: a single
locus held 60.9% of pooled aligned reads in one assembly. Second, the identity was verifiable:
the sequence matched the *P. ostreatus* mitochondrion at 100% identity, and barrnap's
fragmentary rRNA calls fell inside the block. We suggest that clustering unassigned or
intergenic reads and inspecting the top loci should be a routine step when assignment rates are
low, and note that it costs minutes.

### 4.2 The intuitive explanation was wrong, and testing it was cheap

Low gene assignment in 3'-tag data has an obvious candidate explanation: annotated 3' UTRs are
too short, so tag reads fall beyond the last exon. The candidate was well supported here — of
four annotations, only one carried meaningful 3' UTRs, and the others recorded three
nucleotides beyond the stop codon, which is to say none.

Acting on that hypothesis would have been reasonable. Testing it took an afternoon and showed
it to be almost entirely wrong: extending 91.3% of genes by a median of 323 nucleotides moved
assignment from 15.5% to 15.6%. Had we adopted the extension without measuring its effect, we
would have concluded that the annotation was now adequate, retained a plausible-sounding
correction with no benefit, and never looked for the real cause.

The generalisable point is not that UTR extension is useless — it remains correct practice for
tag data, and we retain it — but that the acceptance criterion should be specified before the
correction is applied. The distinction between the two failure modes is visible in the spatial
distribution of unassigned reads, and each implies a different fix.

### 4.3 Reference choice is consequential, and standard metrics do not predict it

The RefSeq-designated reference genome for *P. ostreatus* produced the worst mapping of the
four assemblies tested, by a factor of more than four. Contiguity did not predict performance
either: a contig-level assembly outperformed a chromosome-level one on gene assignment.

The property that mattered is not reported by any standard assembly metric. Assemblies differ
in whether the nuclear rDNA array is resolved into distinct tandem copies or collapsed. Where
it is resolved, rDNA reads accumulate in an identifiable block that can be annotated and
accounted for. Where it is collapsed, the same reads distribute as multi-mappers and are
discarded, inflating apparent multi-mapping and depressing assignment. For expression work on
species with several published assemblies, we therefore recommend empirical selection on
assignment rate using matched annotations, rather than selection on assembly designation,
contiguity or gene count.

That the best-performing assemblies were single-spore isolates of the same cultivar class as
the material studied is unsurprising in hindsight, and reinforces the same recommendation:
strain proximity is worth testing directly rather than assuming that the species reference is
the appropriate reference.

### 4.4 What poly(A) selection did not do

Poly(A) selection is expected to deplete ribosomal RNA, which is not polyadenylated. Here it
did not, and the residual ribosomal fraction varied systematically with tissue. Two aspects are
worth separating. The first is that fungal tissue with high ribosome content presents a harder
depletion problem than the cell types for which such protocols are typically validated; where
sample types of this kind are being profiled, explicit ribosomal depletion may be the safer
choice, and a pilot library is a cheap way to establish which is needed. The second is that
because ribosomal content covaried with tissue identity, it was confounded with the biological
factor of interest, which makes the rRNA-removed analysis track a requirement rather than a
refinement. We could only implement that track correctly once the missing rRNA features had
been added — the two problems compound.

### 4.5 Homology-based annotation systematically misses the interesting proteins

Forty-one per cent of the proteome carried a Swiss-Prot match. That figure sets a hard ceiling
on every downstream enrichment analysis and explains why the two largest co-expression modules
had no enriched terms at all: they are dominated by proteins with no characterised homologue.

For this biological question the bias is unfortunate in a specific way. Lineage-specific
secreted proteins are precisely the class most likely to be involved in a novel exudate-producing
structure, and precisely the class least likely to have a Swiss-Prot homologue. Our secretome
estimate of 576 proteins, obtained by transferring curated signal-peptide annotations from
matched homologues, is a lower bound for the same reason. Absence of enrichment in these data
is therefore weak evidence, and should not be read as absence of function.

### 4.6 What compromised data can still support

Roughly 2.7% of sequenced reads reached a protein-coding feature. It would be defensible to
conclude that such a dataset supports nothing. That conclusion would be wrong, but establishing
where the limit falls required explicit tests rather than judgement.

Two were informative. The first concerns library retention: including the four lowest-yield
libraries made sequencing depth the dominant axis of variation, with the first principal
component correlating with depth at r = 0.686. Removing them eliminated the artefact entirely
(r = −0.003) and *increased* the number of genes passing expression filtering and the number
with a detectable tissue effect. Noise-dominated libraries do not merely add nothing; they
distort dispersion estimation for every other sample. Retaining all samples is not the
conservative choice.

The second concerns which conclusions survive. We applied two criteria sensitive to different
failure modes: a marker test requiring every replicate of a tissue to exceed all samples of
every other tissue, which guards against single-replicate artefacts in a mean-based statistic;
and FDR-corrected co-expression module–trait association, which guards against the multiplicity
inherent in testing dozens of modules. The two agreed completely on which tissues carry signal
and which do not. Where such criteria disagree, neither should be trusted; where they converge,
as here, the conclusion is considerably stronger than either alone.

Both tissues that failed had four and three replicates respectively — more than one that
passed. Sequencing depth, not replicate number, determined what was recoverable.

### 4.7 The metabolic reconstruction

The draft model is offered as a community starting point rather than a predictive tool. Its
construction exposed two limitations worth recording. Mapping enzyme commission numbers onto a
reaction database recovers core metabolism well but systematically misses steps annotated
without a complete EC, which fell disproportionately on cofactor biosynthesis; and transport
reactions, which carry no EC at all, were absent entirely from the initial draft and accounted
for most of its blocked reactions.

Targeted gapfilling resolved the cofactor gaps with three added reactions. One case is
instructive: biotin remained unreachable after gapfilling, and we added it to the medium rather
than to the network, because many fungi are biotin auxotrophs and the culture medium supplies
B-vitamins. Automated gapfilling will readily invent a biosynthetic route for any metabolite
declared essential; whether it should is a biological question, not an algorithmic one.

The model requires curation of mass balance, a biomass composition determined for this organism,
and compartmentalisation beyond cytosol and extracellular space before flux predictions are
meaningful. Context-specific extraction to individual tissues, which was the original intent,
requires expression data of a quality this dataset does not provide.

### 4.8 Limitations

Beyond those already stated: differential expression rests on 12 libraries with two to four
replicates per tissue, and only three libraries exceeded 100,000 assigned counts. The
co-expression network was constructed at a sample size below that recommended by the method's
authors, and its modules are exploratory; we report them with FDR correction and with a
confound test against ribosomal content, but they warrant replication. CAZy assignments are
Pfam-derived class-level calls, not dbCAN family assignments. The secretome is
homology-transferred, not predicted de novo. Conclusions regarding two of the four tissues are
absent not because those tissues lack distinct biology but because these libraries could not
resolve it.

## 5. Conclusions

Three findings generalise beyond this dataset.

First, ribosomal RNA that is absent from a genome annotation is invisible to the rRNA-removal
steps that current consensus pipelines implement, and its loss is reported in a category that
is rarely examined. For *P. ostreatus*, mitochondrial rRNA is unannotated in every public
annotation we examined while constituting the largest single component of these libraries. We
provide corrected annotations and a diagnostic — clustering unassigned reads and asking whether
the loss is focal or diffuse — that distinguishes this from the truncated-UTR explanation that
low assignment in tag data usually attracts.

Second, reference selection for expression analysis should be empirical. The species reference
genome performed worst of four candidates here, and the property that determined performance —
whether the nuclear rDNA array is resolved or collapsed — is not reported by standard assembly
metrics.

Third, severely compromised data can still support conclusions, but only those that survive
criteria sensitive to different failure modes, and only after noise-dominated libraries are
removed rather than retained out of caution. Here two independent criteria converged on the
same two tissues, one of them a newly described structure whose transcriptional signature is
consistently oxidative and secretory.

We provide an rRNA-complete and 3'-aware annotation for two *P. ostreatus* assemblies, a
functional annotation of the proteome including CAZymes and a predicted secretome, and the
first draft genome-scale metabolic reconstruction for the species, achieving non-zero biomass
flux on a defined medium. The complete analysis, including every deviation from the source
pipeline and every correction made during the work, is openly available.

## Data availability
Repository: `https://github.com/dr-richard-barker/Myco_tissue_RNAseq`.
**[TO ADD]** SRA/ENA accession for raw reads; Zenodo DOI for the archived release.

## Author contributions
**[TO WRITE — authors]** CRediT statement.

## Funding
**[TO WRITE — authors]** Grant numbers; must match the `grant_information` field of the
BioProject record in `submission/bioproject.tsv`.

## Acknowledgements
**[TO WRITE — authors]** Consider acknowledging the maintainers of the resources this work
depends on: NASA GeneLab (pipeline specification), Academia Sinica (the *P. ostreatus*
assemblies, including the blue-oyster single-spore isolates that proved the best reference),
UniProt, Pfam, KEGG and ModelSEED.

## Conflicts of interest
**[TO WRITE — authors]** Note that this manuscript reports quality problems in a commercially
supplied dataset. A plain statement of the commercial relationship, or its absence, is
advisable.

## References [TO COMPILE]

Tools and resources requiring citation, grouped by where they appear:

**Culture system and biology** — Porterfield et al. (2026) *Biotechnology Journal*
doi:10.1002/biot.70184.

**Pipeline basis** — NASA GeneLab RNAseq Consensus Pipeline GL-DPPD-7101-G.

**Read processing** — fastp (Chen et al.); HISAT2 (Kim et al.); SAMtools (Danecek et al.);
UMI-tools (Smith et al.); featureCounts / Subread (Liao et al.).

**Statistics** — DESeq2 (Love et al.); edgeR (Robinson et al., Chen et al.); WGCNA (Langfelder
& Horvath); Benjamini & Hochberg (1995); τ specificity index (Yanai et al. 2005).

**Annotation** — barrnap (Seemann); DIAMOND (Buchfink et al.); UniProt Consortium; Pfam
(Mistry et al.); HMMER (Eddy); KEGG (Kanehisa & Goto).

**Metabolic modelling** — ModelSEED (Henry et al., Seaver et al.); COBRApy (Ebrahim et al.).

**Genome assemblies** — Academia Sinica submissions GCA_056149245.1 (BOM_ss5),
GCA_056149315.1 (BOM_ss14), GCA_029852705.2 (PC9.15), GCF_014466165.1 (PC9).

**Suggested context citations to add** — a recent review of white-rot lignocellulose enzymology;
a 3'-tag RNA-seq methods reference (e.g. QuantSeq) establishing the no-length-normalisation
point; a reference for fungal fruiting-body development and cerato-platanin function; and a
reference for ostreolysin/aegerolysin biology in *Pleurotus*.
