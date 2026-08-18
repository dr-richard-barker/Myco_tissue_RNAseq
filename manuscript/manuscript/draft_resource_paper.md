# Ribosomal RNA dominance and unannotated mitochondrial rDNA confound 3′-tag RNA-seq in *Pleurotus ostreatus*: a corrected reference, annotation and metabolic reconstruction

**Draft — resource/methods paper. Target: fungal-biology venue** (e.g. *Fungal Genetics and
Biology*, *Fungal Biology and Biotechnology*, or *IMA Fungus*).

> **Status of this draft.** All numbers are generated from the analysis in this repository and
> are traceable to the files named in each section. Nothing here is estimated or carried over
> from planning documents. Sections marked **[TO WRITE]** need author input. The tissue
> biology is deliberately limited to a demonstration of pipeline function; the full
> tissue-comparison analysis is reserved for a companion paper.

---

## Abstract [TO WRITE — draft below]

Transcriptomic studies of non-model basidiomycetes increasingly rely on 3′-end tag RNA-seq and
on public genome assemblies whose annotations were generated for gene-finding rather than for
read quantification. We show that both choices carry substantial, largely invisible costs. In
16 3′-tag libraries from *Pleurotus ostreatus* cv. "Harbor Blue P01" grown in a mycoponic
ceramic-tube system, ribosomal RNA accounted for the large majority of sequenced material, and
— critically — **the single largest component, mitochondrial rRNA, is unannotated in every
public *P. ostreatus* genome annotation we examined**. Because it carries no feature, it is
neither counted nor removable by biotype filtering, and is silently discarded as "no feature".
A second component, the nuclear rDNA array, is resolved in some assemblies and collapsed in
others, changing multi-mapping rates and therefore quantification. We provide a corrected,
rRNA-complete and 3′-aware annotation, an empirically selected reference, a functional
annotation of the proteome including CAZymes and a predicted secretome, and the first draft
genome-scale metabolic reconstruction for *P. ostreatus*. We further show that a widely
assumed remedy — extending gene 3′ ends to capture tag reads — recovers almost nothing when
the true cause is unannotated rRNA, and we give the diagnostic that distinguishes the two.

---

## 1. Introduction [TO WRITE]

Points to cover:
- *Pleurotus ostreatus* as a white-rot model and a cultivated crop; relevance to controlled-
  environment and space-agriculture bioproduction (cite Porterfield et al. 2026,
  doi:10.1002/biot.70184).
- 3′-end tag counting: cheaper per sample, no length normalisation, increasingly used for
  multi-condition designs. Consequences that are underappreciated — no TPM/FPKM, dependence on
  accurate 3′ ends, sensitivity to polyA-priming artefacts.
- Public fungal annotations are optimised for protein-coding gene models. rDNA is routinely
  omitted; mitochondrial rRNA especially so. State the gap this paper addresses.
- Aim: quantify the cost, provide corrected resources, and give practical diagnostics.

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

### 2.2 Sequencing
Libraries were prepared and sequenced commercially as polyA-selected 3′-end tag libraries
(Illumina NovaSeq X Plus, single-end 94 bp, 14 nt unique molecular identifier appended to the
read name). All 16 libraries were sequenced on one flowcell and lane
(`LH01080:180:25G52WLT4:4`), so no batch term was required. Total yield was 92,430,743 raw
reads (median 5.4 M per library).

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

**This model is a scaffold, not a validated predictor.** Six cofactors (CoA, NAD, NADP, FAD,
biotin, and the chitin precursor) remain unreachable under both references — a systematic
limitation of EC-driven reconstruction, since cofactor biosynthesis steps are frequently
annotated without complete EC numbers. There is no curated biomass objective and no
mass-balance curation. **[PENDING]** Targeted gapfilling of these pathways is in progress.

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

## 4. Discussion [TO WRITE]

Points to develop:
- The mitochondrial rRNA annotation gap is likely general. Fungal genome annotations are
  produced for gene finding; organellar rRNA is routinely absent. Any biotype-based rRNA
  filtering step — including the one in the current GeneLab consensus pipeline — cannot remove
  what is not annotated. Recommend that annotation completeness for rDNA be checked before
  quantification, and give the intergenic-clustering diagnostic (Section 3.4).
- Assembly choice for RNA-seq should weigh rDNA representation, which no standard assembly
  metric reports. The RefSeq reference genome performed worst here.
- 3′-tag protocols: polyA selection did not deplete rRNA in these libraries. Discuss whether
  rRNA depletion should be preferred for basidiomycete tissue with high ribosome content, and
  what QC would have caught this at the provider stage.
- Limits of homology-based functional annotation in basidiomycetes: 41% Swiss-Prot coverage
  means lineage-specific secreted proteins — precisely the interesting ones for exudate
  biology — are systematically invisible.
- The draft GEM as a community starting point, and what curation it needs.

## 5. Conclusions [TO WRITE]

## Data availability
Repository: `https://github.com/dr-richard-barker/Myco_tissue_RNAseq`.
**[TO ADD]** SRA/ENA accession for raw reads; Zenodo DOI for the archived release.

## Author contributions / Funding / Acknowledgements [TO WRITE]

## References [TO COMPILE]
Key citations already required: Porterfield et al. 2026 (doi:10.1002/biot.70184); GeneLab
GL-DPPD-7101-G; HISAT2; fastp; featureCounts/subread; umi_tools; DESeq2; edgeR; WGCNA;
barrnap; DIAMOND; UniProt; Pfam; KEGG; ModelSEED; COBRApy.
