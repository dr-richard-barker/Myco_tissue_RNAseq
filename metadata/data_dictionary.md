# Data dictionary

## `counts/counts_<REF>_{dedup,raw}.txt`
featureCounts output. `dedup` = after UMI deduplication (primary alignments only, umi_tools
directional); `raw` = primary alignments without deduplication. Columns 1–6 are featureCounts
annotation fields; columns 7+ are per-sample counts. Counted with `-s 1` (forward-stranded),
`-t exon -g gene_id`, against `<REF>_final.gtf`.

## `results/read_budget_<REF>.csv`
| column | meaning |
|---|---|
| sample / well / tissue | sample identifiers and factor level |
| raw | reads before trimming (fastp `before_filtering`) |
| rRNA | counts assigned to `gene_biotype "rRNA"` features |
| mRNA | counts assigned to all non-rRNA, non-tRNA features |
| pct | mRNA as a percentage of raw reads |
| det1 / det10 | protein-coding genes with >=1 / >=10 counts |

## `results/dge_<REF>[_filtered]/`
`_filtered` = sensitivity run excluding samples below 15,000 non-rRNA counts (4 dropped).
Each directory holds both an `all_genes` and an `rRNArm` track, independently normalised.
- `normalized_counts_*.csv` — DESeq2 median-of-ratios (linear scale)
- `VST_counts_*.csv` — variance-stabilising transform (log scale; do NOT use for tau)
- `PCA_*.csv` — sample coordinates plus non-rRNA depth, for the depth-artefact check
- `LRT_*.csv` — likelihood-ratio test for any tissue effect
- `DE_<track>_<A>_vs_<B>.csv` — DESeq2 Wald results per pairwise contrast
- `sample_yield.csv` — per-sample non-rRNA counts and low-yield flag

## `results/tissue_models_<REF>/`
- `tau_specificity.csv` — Yanai tau per gene (computed on LINEAR normalised counts) plus
  per-tissue means and the tissue of maximum expression
- `markers_<tissue>.csv` — top 50 genes by tau for that tissue, joined to annotation
- `markers_robust.csv` — the subset where EVERY replicate exceeds the maximum of all samples
  in every other tissue. **Use this set, not the raw marker lists.**

## `results/wgcna/`
- `soft_threshold.csv` — scale-free topology fit across candidate powers
- `gene_modules.csv` — gene to module-colour assignment
- `module_trait_correlation.csv` / `module_trait_pvalue.csv` — Pearson r and p per
  module-eigengene x tissue. **Apply FDR correction: 168 tests.**
- `hub_genes.csv` — top 20 genes per module by module membership

## `results/annotation/`
- `<ref>_vs_sprot.tsv` — raw DIAMOND blastp hits against Swiss-Prot
- `<ref>_functional.tsv` — best hit per protein with EC, GO, KEGG cross-references
- `<ref>_secretome.tsv` — signal peptide and TM features transferred from the Swiss-Prot
  homologue; `secreted` = 1 means signal peptide present and no TM helix beyond it
- `<ref>_pfam.tbl` — HMMER tblout against Pfam-A
- `<ref>_cazymes.tsv` — Pfam domains mapped to CAZy-like functional classes

## `models/`
- `<REF>_draft.xml` — EC-mapped reactions only; GPR rules trace to supporting proteins
- `<REF>_medium.xml` — plus ModelSEED transporters, exchanges constrained to MNM v3
- `<REF>_gapfilled.xml` — plus cofactor gapfilling and a coarse biomass objective.
  Gapfilled reactions have **no gene association** and carry `notes['gapfilled'] = "true"`.

## Identifier note
Count matrices are keyed on GTF `gene_id` (locus tags, e.g. `AAD021_000001`). Annotation
files are keyed on protein accessions (e.g. `KAN1670664.1`). The mapping comes from the
`protein_id` attribute of CDS records in the stock GTF.
