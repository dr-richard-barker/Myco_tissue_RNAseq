# Analysis notes — *Pleurotus ostreatus* tissue-specific RNA-seq (GPNJ7M)

Every deviation from the plan and from NASA GeneLab GL-DPPD-7101-G, with the evidence for it.
Written as the work happened; findings that overturned an earlier assumption are marked.

---

## 1. The vendor analysis is void — wrong species

Plasmidsaurus aligned these reads to *Saccharomyces cerevisiae*.

- `GPNJ7M-expression-matrix.tsv` gene IDs are SGD systematic names: `YDL151C`/BUD30,
  `YER174C`/GRX4, `YHR155W`/LAM1, plus the `Q0045`/COX1, `Q0105`/COB mitochondrial set.
- 7,127 features / 6,600 protein-coding = exactly SGD R64.
- Uniquely mapped: 453–2,284 reads per sample (0.006–0.03% of input).
- Their GSEA used MSigDB Hallmark, which is human-only.

Independent check that the reads themselves are fine: the most abundant read BLASTs at 100%
identity to *Pleurotus* 28S LSU rRNA. Genus confirmed; that region cannot separate
*ostreatus* from *pulmonarius*, so the reference was chosen empirically instead.

Every delivered file — expression matrix, PCA, correlation heatmap, biotype plot, DGE — is
unusable.

## 2. Library chemistry

From the [Plasmidsaurus RNA FAQ](https://plasmidsaurus.com/faq/rna) and the FASTQs:

| Property | Value |
|---|---|
| Assay | 3'-end tag counting (QuantSeq-like), polyA-selected — **not** full-length RNA-seq |
| Reads | SE 94 bp, NovaSeq X Plus, 14 nt UMI appended to read name after `_` |
| Batch | all 16 on one flowcell/run/lane (`LH01080:180:25G52WLT4:4`) — no batch term needed |
| Depth | 2.28M–12.46M raw (median 5.4M) = **11–62% of the advertised 20M** |
| Strandedness | forward (`-s 1`): 27.7% assigned vs 3.1% reverse |
| Pre-processing | vendor already adapter-trimmed (0.002% residual); polyG and polyA remain |

3'-tag counting means **no TPM/FPKM** — one count per molecule regardless of transcript
length, so effective-length normalisation is meaningless. Counts and CPM only.

## 3. STAR is non-functional on this machine — HISAT2 substituted

GL-DPPD-7101-G specifies STAR 2.7.11b. It cannot be used here.

Both the bioconda `osx-64` build (under Rosetta 2) and the native `osx-arm64` build load the
genome correctly (`Loading SA ... done! ... loaded 291869007 bytes`) and then report
`Number of input reads | 0` with `Thread #1 end of input stream, nextChar=-1`.

Ruled out: thread count (fails at `--runThreadN` 1 and 2), `--readFilesCommand`, the sandbox
(fails with it disabled), file location, and the data itself — a synthetic 200 kb genome with
200 perfect reads fails identically. `seqkit` in the same environment reads the same files
without trouble, so it is STAR, not the environment.

HISAT2 2.2.3 passes the same synthetic control at 100% alignment. It is splice-aware and
soft-clips by default. **Deviation accepted and documented.**

## 4. Reference genome: PC9.15

Four candidates, all 16 samples subsampled to 250k reads, genome-only indices:

| Reference | unique% | multi% | unaligned% |
|---|---|---|---|
| BOM_ss5 (`GCA_056149245.1`) | 24.94 | 53.17 | 21.89 |
| BOM_ss14 (`GCA_056149315.1`) | 24.92 | 53.18 | 21.90 |
| PC9.15 (`GCA_029852705.2`) | 21.49 | 53.10 | 25.42 |
| PC9 (`GCF_014466165.1`) | 5.45 | 69.14 | 25.41 |

BOM_ss5 and BOM_ss14 differ by 0.02 points — consistent with being the two nuclei of one
blue-oyster dikaryon. If the culture is that dikaryon, that pair deserves a second look.

**PC9.15 was chosen** on the metric that actually governs counting — gene assignment rate,
12.0% (BOM_ss5) vs 15.5% (PC9.15) on rDNA-depleted reads. It is also the only candidate with
real UTRs and the only chromosome-level assembly:

| Reference | mean 3'UTR | transcripts | genes |
|---|---|---|---|
| PC9.15 | median 78 bp, mean 170 bp | 15,046 | 13,556 |
| PC9 / BOM_ss5 / BOM_ss14 | **3 bp** (stop codon only) | — | 11,849 / 12,705 / 13,310 |

The plan's ">60% unique" gate was written assuming ordinary RNA-seq and does not apply: 53%
of reads are rRNA multi-mapping across repeated rDNA copies. After rRNA depletion, multi-
mapping collapses to 0.4% and unique mapping rises to 49.9%.

## 5. The libraries are dominated by rRNA — from **two** separate sources

The single most important finding for interpreting these data.

**Nuclear rDNA** — measured against barrnap loci, not estimated: **53.5% mean, 19.9–69.7%
range**. The earlier probe-based figure (12–46%, four diagnostic 25-mers) was always a floor
and understated it substantially.

**Mitochondrial rRNA** — **60.9% of everything surviving nuclear depletion**, in one locus:
`CM057219.1:52,604–56,069`, BLAST-confirmed at 100% identity to the *P. ostreatus*
mitochondrion. It is unannotated in all four candidate GTFs (only flanking tRNAs are marked),
so it was invisible until the intergenic reads were clustered. Without adding it, the largest
single component of these libraries cannot be counted at all.

The rRNA share **varies systematically with tissue** (Exudophore 61.7% vs Nodule 45.5%), so
it is confounded with the factor of interest. This is why the GL-DPPD-7101-G rRNA-removed
track is mandatory here rather than optional.

## 6. Corrected assumption: the 3'UTR extension was *not* the problem

The plan predicted that missing 3'UTRs were causing the 30–78% "no feature" rate, and that
extending gene 3' ends would recover it. It was built and measured:

- 12,380 of 13,556 genes extended (91.3%), median 323 bp, neighbour-capped at 50 bp.
- Effect on assignment: **15.5% → 15.6%**.

The plan's own acceptance check said to treat a non-material rise as a signal that the logic
was wrong. It was — the real cause was the unannotated mitochondrial rRNA above. The
extension is retained (it costs nothing and is correct for a 3'-tag assay) but it is not the
fix that was predicted.

## 7. Pipeline deviations from GL-DPPD-7101-G

| GeneLab | Here | Reason |
|---|---|---|
| GL-DPPD-7101-F (`NF_RCP-F`) | GL-DPPD-7101-**G** | G is current and adds the rRNA-removed DGE track |
| STAR 2.7.11b | HISAT2 2.2.3 | STAR non-functional on this OS (§3) |
| RSEM | featureCounts (subread 2.1.1) | RSEM models full-length coverage; invalid for 3'-tag |
| TrimGalore! | fastp 1.3.6 | reads pre-trimmed; the real issues are polyG and polyA |
| Ensembl auto-fetch | NCBI, chosen empirically | *P. ostreatus* is not in Ensembl |
| annotation as-is | PC9.15 + 3' extension + nuclear/mito rRNA features | §5, §6 |
| TPM/FPKM reported | counts/CPM only | no effective-length normalisation for 3'-tag |
| paired-end | single-end | SE 94 bp library |
| Nextflow + Singularity | conda envs, versions pinned to G | Singularity has no native macOS support |
| — | umi_tools directional dedup | 14 nt UMI present and essential for 3'-tag |

DESeq2 1.46.0, samtools 1.21, subread 2.1.1, FastQC 0.12.1 all match the G spec exactly.

## 8. Implementation gotchas worth remembering

- **fastp picks compression from the filename.** Writing to `*.fq.gz.part` produced a plain
  file that was then renamed to `.gz`; HISAT2 silently aligned 0 reads across all 16 samples.
  Temp names must keep the `.gz` suffix, and the pipeline now runs `gzip -t` before renaming.
- **Annotating a locus on both strands makes reads ambiguous under `-s 0`.** The first
  mitochondrial rRNA annotation used `+` and `-`, which silently discarded 27k–102k reads per
  sample as `Unassigned_Ambiguity`. One feature per region now.
- **Secondary alignments crushed umi_tools.** HISAT2 emitted 8.9M secondary records per 2.5M
  primary (rRNA multi-mapping). Dedup projected to 3.7 h; filtering to primary alignments
  first (`-F 0x104`, which featureCounts discards anyway) cut it to ~20 min.
- **`micromamba run` per command serialises on a package-cache lock.** Envs go on `PATH`.
- **zsh does not word-split unquoted variables** — `$MM` holding a multi-word command fails.
- **`set -o pipefail` + `head` closing a pipe** kills `gzcat` with SIGPIPE and aborts the script.

## 9. Functional annotation: Swiss-Prot instead of eggNOG

`eggnogdb.embl.de` does not resolve from this machine (UniProt, EBI, NCBI and GenomeNet all
do). Substituted DIAMOND against Swiss-Prot:

- 6,158 of 14,901 PC9.15 proteins hit (41%)
- 3,439 proteins carry at least one EC number
- 1,929 distinct EC numbers

## 10. Draft GEM status

`models/PC9.15_draft.xml` — 2,055 reactions, 2,222 metabolites, 1,042 genes, built by
EC-mapping the annotated proteome onto the ModelSEED biochemistry database. GPR rules trace
every reaction to the supporting PC9.15 proteins.

**This is a raw scaffold, not a usable model.** Measured quality:

| Metric | Value |
|---|---|
| blocked reactions (cannot carry flux) | 1,688 / 2,071 (81.5%) |
| reactions that can carry flux | 383 |
| mass-unbalanced reactions | 515 (24.9%) |
| metabolites with no formula | 140 (6.3%) |
| metabolites in exactly one reaction | 1,221 (55.0%) |
| exchange reactions | 16 |

Outstanding before any flux prediction is meaningful: transport reactions (ModelSEED
transporters carry no EC so none were captured, which is the main cause of the 81.5% blocked
rate), a biomass objective, gapfilling against the real culture medium (**still needed from
the user**), compartmentalisation beyond cytosol/extracellular, and mass-balance curation.

### 10b. GEM on the real culture medium

The culture medium was supplied by the user: mycoponic nutrient medium (MNM v3) from
Porterfield et al. 2026 ([10.1002/biot.70184](https://doi.org/10.1002/biot.70184)) — see
`metadata/culture_conditions.md`. Adding ModelSEED transport reactions and constraining
exchanges to that medium (`models/PC9.15_medium.xml`):

| Metric | Draft | On MNM |
|---|---|---|
| reactions | 2,071 | 5,409 |
| able to carry flux | 383 | **2,280** |
| blocked | 81.5% | 57.8% |
| biomass precursors synthesisable | — | **28 / 34** |

All 20 proteinogenic amino acids and all NTPs/dNTPs are reachable from the medium. **See §14
for a correction to the cofactor result reported here** — two of the six ModelSEED compound
ids used for that test were wrong.

This is a connectivity diagnostic, not a validated growth prediction: there is still no
curated biomass objective and no mass-balance curation.

Per-tissue contextualisation (RIPTiDe/Troppo) remains deferred until the counts are trusted
and the cofactor gaps are filled — extracting context-specific models from a scaffold that
cannot make NAD would produce confident-looking nonsense.

## 11. Results

### 11a. Read budget (full depth, all 16 samples)

Only **3 of 16** samples exceed 100k assigned non-rRNA counts: 2H (1,250,686), 1A (448,440),
1G (204,742). The rest span 2,454–74,228. The dominant loss is multi-mapping (45.6–94.7% per
sample) — nuclear rDNA hitting repeated copies, which featureCounts discards. Strandedness
re-verified on full-depth BAMs: `-s 1` 49.6%, `-s 0` 47.6%, `-s 2` 6.3%.

Full table in `results/read_budget.csv`.

### 11b. Four samples were poisoning the analysis

The PCA showed Fuzzy mycelium clustering tightly across a 5.5x depth spread, and Nodule
(2H/2F/2G) clustering tightly across a **31x** spread — real biology. The scattered samples
were exactly the four lowest-yield ones: 1C (1,594), 2E (4,410), 2B (6,695), 2A (8,107).

Dropping those four (`results/dge_filtered/`, threshold 15k):

| | all 16 | 12 retained |
|---|---|---|
| cor(PC1, log10 depth) | 0.686 | **0.014** |
| genes passing filterByExpr (rRNArm) | 1,079 | 1,541 |
| LRT any-tissue effect (padj<0.05) | 76 | **261** |
| Nodule vs Fuzzy mycelium | 81 | **227** |
| Nodule vs Exudophore | 70 | **168** |

The depth artefact disappears completely. All four tissues retain n>=2 (Exuding 3,
Exudophore 2, Fuzzy 4, Nodule 3). Both runs are kept: `results/dge/` (all 16, per the
user's instruction to flag rather than drop) and `results/dge_filtered/` (sensitivity).

`Fuzzy mycelium vs Exuding mycelium` returns **0 DE genes in both runs** — consistent, and
plausibly real given both are vegetative mycelium.

### 11c. Per-tissue markers are biologically coherent

70 genes at tau >= 0.85 (Nodule 47, Exudophore 12, Fuzzy 9, Exuding 2). Top annotated markers:

- **Nodule** — small heat shock proteins (17.6/17.8 kDa class I, HSP16) and *cerato-platanin*
  family ("eliciting plant response-like protein CP1"). Consistent with nodules as early
  fruiting-body primordia.
- **Exudophore** — an oxidative/secretory signature: **glyoxal oxidase (GLOX)**, copper amine
  oxidase, formate and aldehyde dehydrogenases, peroxisomal acyl-activating enzyme, secreted
  zinc metalloprotease. GLOX and copper amine oxidase both generate extracellular H2O2, the
  co-substrate for white-rot peroxidases.
- **Fuzzy mycelium** — **ostreolysin A6** (an aegerolysin named for *P. ostreatus*, linked to
  fruiting initiation), fruiting body differentiation protein 16, chitinase 1, beta-1,3-glucan
  binding protein, and an **LPMO** (AA9 lignocellulose CAZyme).
- **Exuding mycelium** — aryl-alcohol dehydrogenase and GMC oxidoreductase (lignin-related
  white-rot machinery), O-methyltransferase, ammonium transporter 1.

The LPMO and aryl-alcohol dehydrogenase signals fit the oak-sawdust/cellulose medium
(§ culture conditions). **Caveat:** only 1,541 genes survived filtering, so this is a narrow
window, and tau rests on 4 tissue means with n=2-4.

**Bug caught:** tau was first computed on VST values and returned 0 tissue-specific genes —
log compression drives every value/max ratio to 1. It must be computed on linear normalised
counts. Separately, the marker/annotation join silently returned nothing because counts are
keyed on GTF `gene_id` (`PTI98_*`) while the annotation is keyed on protein accessions
(`KAJ*`); the map comes from the GTF CDS records (13,411 genes, 5,385 annotated).

### 11d. Reference re-test on matched annotations

Promised in §10c. Both annotations 3'-extended, same trimmed reads, same settings, 4 samples:

| Reference | total assigned | multi-mapping (range) |
|---|---|---|
| **BOM_ss5** | **5,667,553** | 23.1–79.4% |
| PC9.15 | 5,331,061 | 25.1–81.3% |

BOM_ss5 assigns **6.3% more reads** and has consistently lower multi-mapping — as expected
now the cultivar is confirmed blue oyster. The margin is modest and would not change any
conclusion above; PC9.15 remains structurally the better assembly (13 chromosomes vs 16
contigs, 13,556 genes with isoforms vs 12,705). Re-running on BOM_ss5 is scripted and cheap
(~1 h) if the strain match is judged more important than assembly contiguity.

## 12. Full re-run on BOM_ss5 (final reference)

On the user's instruction the entire analysis was re-run against `BOM_ss5`
(`GCA_056149245.1`), the reference matching the confirmed cultivar. Trimmed FASTQs are
reference-independent and were reused; `10_pipeline.sh` now takes a reference label and
writes to `bam_<label>/` and `counts/counts_<label>_*.txt`.

### 12a. BOM_ss5 resolves the rDNA array that PC9.15 collapses

Clustering pooled intergenic reads found a repeating block structure at
`JBQVBD010000012.1:2,298,112-2,346,320` (**40.2%** of pooled aligned reads) plus
mitochondrial rRNA on `CM148777.1` (**31.3%**) — ~71% in identifiable rDNA loci. barrnap's
fragmented 5S/5.8S/18S/28S calls fall inside every block, confirming them. PC9.15 collapses
the nuclear array, so those reads scatter into multi-mapping instead. That is the mechanism
behind BOM_ss5's higher assignment rate.

rRNA blocks are now declared per reference in `refs/rRNA_regions_<label>.tsv` rather than
hardcoded.

### 12b. BOM_ss5 wins on every sample

Assigned non-rRNA counts, **+6.9% overall, and higher for all 16 samples** (range +4.0% to
+13.5%): total 2,316,587 (PC9.15) -> 2,475,808 (BOM_ss5). Best samples: 2H 1,250,686 ->
1,339,504; 1A 448,440 -> 469,981; 1G 204,742 -> 224,213.

The sample-level verdict is unchanged: the same 3 samples clear 100k counts and the same 4
are noise-dominated.

### 12c. DE improves consistently (12 retained, rRNA-removed track)

| | PC9.15 | BOM_ss5 |
|---|---|---|
| genes passing filterByExpr | 1,541 | **1,666** |
| cor(PC1, log10 depth) | 0.014 | **-0.003** |
| LRT any-tissue effect | 261 | 281 |
| Nodule vs Exudophore | 168 | **198** |
| Nodule vs Fuzzy mycelium | 227 | 250 (all-genes track) |
| tissue-specific genes (tau>=0.85) | 70 | **82** |

`Fuzzy vs Exuding mycelium` remains ~0-2 DE, reproducing across both references — good
evidence it is a real biological similarity rather than a reference artefact.

### 12d. The marker biology replicates independently

Running against a different genome and a separately annotated proteome reproduced every
tissue signature, which is the strongest available evidence they are real:

- **Nodule** — cerato-platanin (CP1, twice), small HSPs (HSP16, 17.8 kDa class I), PRY1
  (CAP-superfamily secreted sterol-binding protein). Fruiting-body primordium signature.
- **Exudophore** — an even clearer extracellular-H2O2 signature than under PC9.15:
  **alcohol oxidase 1**, **glyoxal oxidase (GLOX)**, FAD-dependent monooxygenase, formate
  dehydrogenase, acetate--CoA ligase, amino-acid permease. Alcohol oxidase and GLOX both
  generate the H2O2 that white-rot peroxidases require.
- **Fuzzy mycelium** — **ostreolysin A6**, **dye-decolorizing peroxidase (DyP)** — a genuine
  lignin-modifying peroxidase, newly surfaced by this reference — fruiting body
  differentiation protein 16, LPMO (AA9), extracellular metalloprotease.
- **Exuding mycelium** — aryl-alcohol dehydrogenase, O-methyltransferase tpcA, versicolorin
  reductase, ammonium transporter 1, ABC transporter.

### 12e. BOM_ss5 GEM

| | PC9.15 | BOM_ss5 |
|---|---|---|
| proteins | 14,901 | 12,521 |
| Swiss-Prot hits | 6,158 | 5,163 |
| distinct EC numbers | 1,929 | 1,795 |
| draft reactions | 2,055 | 1,899 |
| flux-carrying on MNM | 2,280 | 2,262 |
| biomass precursors | 28/34 | 28/34 |

`models/BOM_ss5_draft.xml`, `models/BOM_ss5_medium.xml`. **Superseded by §14** — the cofactor
comparison reported here used two incorrect compound ids.

PC9.15's GEM is marginally richer because its annotation carries isoforms (14,901 vs 12,521
proteins). For transcript quantification BOM_ss5 is better; for the metabolic scaffold the
two are near-equivalent.

### 10c. Strain identity

The cultivar is ***P. ostreatus* cv. "Harbor Blue P01"**, a commercial **blue oyster**
(Porterfield et al. 2026). This is relevant to §4: `BOM_ss5`/`BOM_ss14` are Academia Sinica
**B**lue **O**yster **M**ushroom single-spore isolates and mapped 3.5 points better than
PC9.15 on raw unique rate. PC9.15 won the assignment comparison only because it was the sole
candidate carrying real UTRs. With the 3' extension now applied to both (BOM_ss5: 97.7% of
genes extended, median 500 bp), that advantage is neutralised and the comparison should be
re-run on matched annotations before the reference choice is considered final.

## 13. Marker robustness — a correction

tau is computed on tissue **means**, so a single outlying replicate can manufacture a
high-tau "marker" for a whole tissue. `scripts/18_marker_robustness.py` applies the test the
mean cannot: a marker counts only if **every** replicate of its tissue exceeds the highest
value in any sample of any other tissue.

| Tissue | n | top-50 markers supported by all replicates |
|---|---|---|
| Exudophore | 2 | 37/50 (74%) |
| Nodule | 3 | 22/50 (44%) |
| Exuding mycelium | 3 | 1/50 (2%) |
| Fuzzy mycelium | 4 | **1/50 (2%)** |

**This corrects §12d.** Re-running against a second reference reproduced the marker lists,
and I described that as the biology "replicating independently". That was wrong: both
references were quantifying *the same reads*, so agreement tests reference robustness, not
biological reproducibility. The across-replicate test above is the real one, and it is much
less kind.

What survives:

- **Nodule** — holds up. Cerato-platanin (CP1) and the small HSPs are supported by 2F, 2G and
  2H independently (HSP16: 220/214/1399 vs 103 max elsewhere), so this is not a 2H artefact
  despite 2H carrying 94% of the group's counts.
- **Exudophore** — alcohol oxidase 1, formate dehydrogenase, acetate--CoA ligase, FAD
  monooxygenase and an AA9 LPMO are consistent across both replicates. But n=2 makes this
  test far easier to pass, so 74% overstates the confidence. **Glyoxal oxidase (GLOX) does
  NOT pass** (118/613 vs 294 elsewhere) — the H2O2-generation story rests on alcohol oxidase,
  not GLOX as previously emphasised.
- **Fuzzy mycelium** — does **not** hold up. Every top marker is driven by one replicate,
  mostly 1F. This includes the **dye-decolorizing peroxidase** highlighted in §12d as a new
  finding, and ostreolysin A6 (98/192/30/236 vs 92 elsewhere — 1G falls below). Treat the
  entire fuzzy-mycelium signature as unsupported.
- **Exuding mycelium** — only 2 tau-specific genes existed; 1/50 survives. Nothing to report.

Robust markers only: `results/tissue_models_BOM_ss5/markers_robust.csv` (61 genes).

Note the inverse relation between n and pass rate (n=2 -> 74%, n=4 -> 2%): the test gets
harder with more replicates, which is the correct behaviour but means cross-tissue
comparison of these percentages is not meaningful.


## 14. Correction: two cofactor ids were wrong, and the gapfilling method was intractable

Two separate problems in the metabolic work, found while redoing the gapfilling.

### 14a. Wrong compound identifiers

The blocked-cofactor test used `cpd00166` for CoA and `cpd00557` for the chitin precursor.
Checked against `refs/modelseed/compounds.tsv`, those ids are **Calomide** (a cobalamin,
C72H100CoN18O17P) and **Siroheme** (C42H36FeN4O16). Neither is CoA or a chitin precursor.

The correct ids are **CoA = `cpd00010`** and **UDP-N-acetylglucosamine = `cpd00037`**. NAD
(`cpd00003`), NADP (`cpd00006`), FAD (`cpd00015`) and biotin (`cpd00104`) were right.

This error was introduced in `14_build_gem.py` / `16_gem_medium.py` and propagated into the
notes, the README and the manuscript draft as "six blocked cofactors: CoA, NAD, NADP, FAD,
biotin, chitin precursor". That statement was wrong: CoA was never blocked, and the two
compounds actually being tested were unrelated to the claim.

### 14b. The global MILP gapfill was the wrong method

The first attempt passed the entire ModelSEED universal pool (~31k candidate reactions) to
COBRApy's `gapfill()` for each target. That makes MILP size a function of the database rather
than the question; it ran **10.5 h without converging on a single target** and was killed.

`scripts/26_gem_gapfill_targeted.py` restricts the candidate pool by walking backward through
the universal network from each target metabolite (depth 6, capped), which reduces the pool to
1.8k-8k reactions. Every target then solves in **2.5-9.3 s**.

### 14c. Corrected result

With the right ids and the tractable method (`models/BOM_ss5_gapfilled.xml`):

| Target | Outcome |
|---|---|
| CoA | already producible — was never blocked |
| NAD | gapfilled, 1 reaction (`rxn07095`) |
| NADP | already producible |
| FAD | gapfilled, 1 reaction (`rxn00122`) |
| UDP-GlcNAc (chitin precursor) | gapfilled, 1 reaction (`rxn28025`) |
| biotin | **supplied by the medium, not synthesised** |

Biotin remained unreachable after targeted gapfilling, and that is biologically correct rather
than a modelling gap: many fungi are biotin auxotrophs, and MNM v3 contains malt extract,
peptone and tryptic soy broth, all of which supply B-vitamins. It was therefore added to the
medium definition rather than given a fabricated biosynthesis route.

Final state, both references, on the MNM medium:

| | BOM_ss5 | PC9.15 |
|---|---|---|
| reactions | 5,247 | 5,414 |
| able to carry flux | 2,287 | 2,305 |
| blocked | 56.4% | 57.4% |
| cofactors producible | **6/6** | **6/6** |
| biomass flux on MNM | **125.0** | **125.0** |
| gapfilled reactions (no gene evidence) | 3 | 3 |

The models now produce non-zero biomass flux, which the earlier versions did not. The biomass
objective remains coarse and uncurated, and gapfilled reactions carry no gene association and
are tagged `notes['gapfilled']='true'`.
