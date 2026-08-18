# Tissue-specific transcriptomics of *Pleurotus ostreatus* mycoponic mycelium

Reanalysis of 16 3'-end tag RNA-seq libraries from four mycelial tissue types of
*Pleurotus ostreatus* cv. "Harbor Blue P01", grown in a mycoponic ceramic-tube system
(Porterfield et al. 2026, [doi:10.1002/biot.70184](https://doi.org/10.1002/biot.70184)).

Pipeline derived from the NASA GeneLab RNAseq Consensus Pipeline
**[GL-DPPD-7101-G](https://github.com/nasa/GeneLab_Data_Processing/blob/master/RNAseq/Pipeline_GL-DPPD-7101_Versions/GL-DPPD-7101-G.md)**,
with every deviation documented and justified in [`NOTES.md`](NOTES.md).

---

## Why this reanalysis exists

The sequencing provider's delivered analysis aligned these reads to
***Saccharomyces cerevisiae***. Gene IDs in the delivered matrix are SGD systematic names
(`YDL151C`, `YER174C`, `Q0045`/COX1); 7,127 features / 6,600 protein-coding is exactly SGD
R64; uniquely mapped reads were 453–2,284 per sample (0.006–0.03% of input). The delivered
expression matrix, PCA, correlation heatmap, biotype plot and DGE results are unusable.

The reads themselves are fine — the most abundant read BLASTs at 100% identity to
*Pleurotus* 28S rRNA.

## Headline findings

1. **The libraries are dominated by rRNA from two separate sources.** Nuclear rDNA is
   19.9–69.7% of reads (mean 53.5%); mitochondrial rRNA accounts for a further ~60% of what
   survives nuclear depletion. Neither is annotated in any public *P. ostreatus* GTF, so both
   had to be located from coverage and added to the annotation.
2. **Only 3 of 16 libraries exceed 100k assigned non-rRNA counts.** Effective mRNA yield is
   0.1–10.5% of raw reads. Four libraries are noise-dominated.
3. **Excluding those four removes a severe depth artefact** — `cor(PC1, depth)` falls from
   0.686 to −0.003 — and roughly triples the number of differentially expressed genes.
4. **Two tissues carry defensible signal: Nodule and Exudophore.** Both a replicate-supported
   marker test and FDR-corrected WGCNA module–trait association independently identify these
   two and only these two. Fuzzy mycelium and Exuding mycelium do not survive either test.

## Experimental design

| Tissue | Wells | Files | n |
|---|---|---|---|
| Exuding mycelium | 1A–1D | `GPNJ7M_1`–`_4` | 4 |
| Fuzzy mycelium | 1E–1H | `GPNJ7M_5`–`_8` | 4 |
| Exudophore | 2A–2D | `GPNJ7M_9`–`_12` | 4 |
| Nodule | 2E–2H | `GPNJ7M_13`–`_16` | 4 |

Single factor, 4 levels. All 16 libraries on one flowcell/lane — no batch term required.
Library chemistry: Illumina NovaSeq X Plus, single-end 94 bp, polyA-selected 3'-end tag
counting, 14 nt UMI in the read name. **3'-tag counting means no TPM/FPKM** — counts and CPM
only, since one count represents one molecule regardless of transcript length.

## Reference genome

Chosen empirically, not assumed. Four candidates were test-mapped; the final analysis uses
**BOM_ss5 (`GCA_056149245.1`)**, a Blue Oyster Mushroom single-spore isolate matching the
cultivar. It assigns 6.9% more reads than the runner-up on every one of the 16 samples, and
uniquely resolves the tandem nuclear rDNA array that other assemblies collapse.

A complete parallel analysis against **PC9.15 (`GCA_029852705.2`)** is retained
(`results/dge_*`, `results/tissue_models/`, `models/PC9.15_*`) so the choice is auditable.

## Repository layout

```
scripts/     numbered, run in order; each documents its own rationale
metadata/    runsheet, culture conditions, reference provenance
counts/      featureCounts matrices (per reference, dedup + raw)
results/     read budgets, DGE, tissue models, WGCNA, annotation
models/      SBML metabolic models (draft, medium-constrained, gapfilled)
envs/locks/  explicit conda lockfiles for every environment
logs/        run logs
NOTES.md     every deviation and correction, with evidence
```

Large regenerable artefacts (BAMs, trimmed FASTQs, reference databases, conda environments;
~37 GB) are excluded via `.gitignore`. They are reproducible from `scripts/` plus
[`metadata/reference_sources.tsv`](metadata/reference_sources.tsv).

## Reproducing

```bash
bash scripts/00_setup_envs.sh          # conda environments
bash scripts/01_fetch_refs.sh          # candidate genomes from NCBI
bash scripts/02_testmap.sh             # empirical reference selection
python3 scripts/03_testmap_report.py
bash scripts/04_rrna_ref.sh BOM_ss5    # barrnap rRNA loci
python3 scripts/07_extend_3p.py refs/BOM_ss5/BOM_ss5_genomic.gtf refs/BOM_ss5/BOM_ss5_3p500.gtf
python3 scripts/08_build_annotation.py refs/BOM_ss5/BOM_ss5_3p500.gtf refs/rRNA/BOM_ss5_rRNA.bed \
        refs/BOM_ss5/BOM_ss5_final.gtf --regions refs/rRNA_regions_BOM_ss5.tsv
bash scripts/10_pipeline.sh BOM_ss5    # trim, align, dedup, count
python3 scripts/11_counts_report.py --counts counts/counts_BOM_ss5_dedup.txt --gtf refs/BOM_ss5/BOM_ss5_final.gtf
Rscript scripts/12_dge.R counts/counts_BOM_ss5_dedup.txt 15000 drop BOM_ss5
python3 scripts/15_tissue_models.py    # tau specificity + markers
python3 scripts/18_marker_robustness.py
Rscript scripts/20_wgcna.R
```

Raw FASTQ files are not included; they are the property of the data generator.

## Known limitations

- **Underpowered.** 12 of 16 libraries retained; tissue n = 2–4. WGCNA is run at n=12, below
  the ~15 its authors recommend — modules are exploratory and reported with FDR correction.
- **Fuzzy mycelium and Exuding mycelium conclusions are not supported** by either robustness
  test and should not be reported as findings.
- **Secretome is homology-transferred, not predicted.** SignalP/Phobius are licence-restricted
  and DeepSig's conda build pins an unavailable TensorFlow; signal peptides are inherited from
  the best-hit Swiss-Prot entry, covering only ~41% of the proteome. Counts are a lower bound.
- **CAZyme assignment is Pfam-derived, not dbCAN.** The dbCAN download service has moved
  behind a landing page and its database files were unreachable.
- **The metabolic model is a draft.** It carries an uncurated coarse biomass objective, three
  gapfilled reactions with no gene evidence (tagged `notes['gapfilled']`), and no mass-balance
  curation. It produces non-zero biomass flux on the defined medium and all six tested
  cofactors are producible, but it is a scaffold for hypothesis generation, not a validated
  flux predictor. An earlier version of these notes reported six blocked cofactors; two of the
  compound ids in that test were wrong — see `NOTES.md` §14.
- **No PAS atlas.** 3' ends were extended by a flat neighbour-capped rule rather than called
  from polyadenylation sites; the extension changed gene assignment by only 0.1%, so the
  planned PAS refinement was not pursued.

## Citation

See [`CITATION.cff`](CITATION.cff).

## Licence

Code: MIT (see [`LICENSE`](LICENSE)). Reference data retain their original licences.
