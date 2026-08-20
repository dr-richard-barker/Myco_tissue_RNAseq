# Documentation audit

A line-by-line check of the six manuscripts, `NOTES.md` and `README.md` for fabricated,
mis-attributed or stale claims. Prompted by a reconnaissance pass that found a fabricated
bibliography entry; widened once it became clear that nothing in the repository checked
citations, identifiers or cross-document consistency, and that the one numeric check that did
exist had blind spots wide enough to hide most of what is listed below.

Everything here has been corrected in the sources unless explicitly marked otherwise.

---

## 1. A fabricated citation

`10.3389/fmicb.2020.00863` — the paper the collaboration supplied as the basis for the
mitochondrial work — was recorded as:

> Zaccaron AZ, Bluhm BH. *Erosion of Genome Sequences and Emergence of Novel Genes in the
> Mitochondrial Genome of the Fungal Pathogen* Bipolaris sorokiniana.

The DOI is correct. Everything attached to it was invented. The real paper is:

> Song N, Geng Y, Li X (2020). *The Mitochondrial Genome of the Phytopathogenic Fungus*
> Bipolaris sorokiniana *and the Utility of Mitochondrial Genome to Infer Phylogeny of
> Dothideomycetes.* Frontiers in Microbiology 11: 863.

Alex Zaccaron is a real author of fungal mitogenome papers, including a 2021 grape powdery
mildew mitogenome. This is therefore a conflation of a real researcher with the wrong paper —
the hardest kind of error to notice, because every component looks plausible in isolation.

Corrected in `latex/refs.bib` (key renamed `song2020bipolaris`) and at all five in-text uses in
`P4_mitochondrial/sections.tex` and `P5_mitocomparative/sections.tex`, including the prose
attribution "Zaccaron and Bluhm's analysis of…".

**Not corrected:** the commit message of `38a93f3` carries the same wrong attribution.
Rewriting the message would rewrite published history on a public repository for no benefit to
the record. It is documented here instead. This is the maintainer's call to revisit.

## 2. Claims made *about* that paper

P4's introduction asserted the source describes "protein-coding genes interrupted by numerous
group I and group II introns, a full tRNA complement, and dozens of free-standing ORFs of
unknown function". The group I / group II attribution is **not** supported: the paper reports
28 introns without classifying them by type.

Replaced with what the paper does report — a 137,775 bp genome, 12 core protein-coding genes,
28 introns, 38 tRNAs and 52 free-standing ORFs of unknown function.

P5's introduction attributed to the same paper the general claim that group I introns and their
homing endonucleases expand fungal mitogenomes. That claim is real but belongs to other work.
Re-attributed to two sources that make it, both verified against CrossRef and PubMed:

- Sandor, Zhang & Xu (2018) *Appl Microbiol Biotechnol* 102:9433–9448 — fungal mitogenome
  structure, size and gene content.
- Férandon, Xu & Barroso (2013) *Fungal Genet Biol* 55:85–91 — the 135 kbp *Agaricus bisporus*
  mitogenome, 43 group I introns occupying 45% of the sequence. An agaric, so a closer
  comparator for *Pleurotus* than the ascomycete was.

The same two now source P4's previously unsourced generalisation about what a typical fungal
mitogenome encodes.

## 3. Bibliography completeness

- `porterfield2026mycoponics` listed one author of ten. Completed from CrossRef. This also
  settled an open question about the manuscript author list: the published record gives
  **Sanchez, Adriana K.** and **Moulton, Simone X.**, confirming those given-name/surname
  splits, and lists the corresponding author of this work as **Richard J. Barker**.
- `benjamini1995controlling` had no DOI. Added.
- `danecek2021samtools` was in the bibliography but never cited; SAMtools does the primary-
  alignment filtering, so it is now cited in the Methods where that step is described.
- All 28 remaining entries verified: title, full author list, year, volume and first page.

## 4. Tool versions stated but never checked

`NOTES.md` claimed "DESeq2 1.46.0, samtools 1.21, subread 2.1.1, FastQC 0.12.1 all match the G
spec exactly." Two of those version numbers were the ones named in GL-DPPD-7101-G, not the ones
installed. The environments contain exactly one DESeq2 (**1.50.2**) and one SAMtools (**1.24**),
so those are what produced the results.

Corrected in `NOTES.md` and in the Methods of P1 and V0. The pipeline steps are equivalent to
the GeneLab spec; the versions are newer, and that is now stated rather than glossed. FastQC is
installed but unused — fastp performs the QC — which is also now stated.

## 5. Stale numbers left behind by the switch to BOM_ss5

Three findings share one cause: results were regenerated on BOM_ss5, but figures computed
during the earlier PC9.15 work were not all revisited.

| Claim | Was | Is | Source of truth |
|---|---|---|---|
| PC1 correlation with depth, all 16 libraries | `r = 0.686` | `r = 0.681` | 0.686 is the PC9.15 value in `results/dge/`; 0.681 is BOM_ss5 |
| Nuclear rDNA array share of pooled aligned reads | 40.2% | **59.7%** | direct recount over the 16 pooled primary BAMs |
| Mitochondrial contig share | 31.3% | **29.1%** | as above; agrees with `results/figure_data/intergenic_peaks.csv` |
| Combined "identifiable rDNA loci" | ~71% | **~89%** | as above |

The corrected rDNA figures are *stronger* than the ones claimed: nearly nine in ten aligned
reads sit in identifiable ribosomal loci, not seven in ten. The recount matches the independent
`mito total alignments` fact (21,180,726) exactly, which is what gives confidence in it.

## 6. An experiment attributed to the wrong reference

The 3′-extension result — "91.3% of 13,556 genes extended; median extension 323 bp; 3.79 Mb
added" — was presented in the Methods of a BOM_ss5 pipeline. Those are PC9.15 numbers:
BOM_ss5 has 12,705 genes, PC9.15 has 13,556, and re-running `scripts/07_extend_3p.py` on each
reproduces 91.3%/323 bp/3.79 Mb for PC9.15 and 97.7%/500 bp/4.55 Mb for BOM_ss5.

The experiment itself was legitimately done on PC9.15 — Supplementary Fig. S2 is built from
`qc/fcmp/PC9.15*`, and PC9.15 is the only annotation with meaningful UTRs, which makes it the
fairest test of the hypothesis. The error was that the text never said so. The passage now
attributes the test to PC9.15, explains why that is the favourable case, and adds the BOM_ss5
equivalent. The figure caption now names PC9.15 too.

## 7. An ambiguous sample set

The rRNA-versus-depth correlation of `r = -0.907` is computed over the twelve libraries
entering the co-expression network, not all sixteen — over all sixteen it is −0.83. The text
said "these libraries" in a context where either reading was available. Now stated explicitly.

## 8. Not reproducible, and marked as such

The claim that a single mitochondrial locus held **60.9%** of pooled aligned reads under PC9.15
cannot be recomputed: the PC9.15 BAMs were not retained after the switch to BOM_ss5. It is
consistent with `NOTES.md` written at the time, and it is reported as a PC9.15 observation, so
it stays — but it rests on a contemporaneous note rather than on a surviving artefact.

## 9. Bugs in the audit tooling itself

`scripts/32_audit_numbers.py` was the only pre-existing check on quantitative claims. It was
substantially less effective than it appeared:

- **Escaped percent signs truncated the text.** Comments were stripped with `%.*`, which also
  matches the `%` in `59.7\%` and deleted the rest of that line. Every claim appearing after
  the first percentage on a line was invisible to the audit — a large fraction of the results
  sections.
- **`\texttt{}` with inner braces leaked digits.** LaTeX digit grouping
  (`\texttt{JBQVBD010000012.1:2{,}298{,}112--2{,}346{,}320}`) terminated the `[^}]*` argument
  match at the first inner brace, so 298, 112, 320 and 346 were audited as if they were claims.
- **Layout was audited as data.** `\includegraphics[width=0.62\textwidth]` and
  `\parbox{0.97\textwidth}` contributed 0.62 and 0.97 as quantities.
- **Version strings became quantities.** "HISAT2 2.2.3" was read as the number 2.2.
- **The number regex could not parse `2{,}135`**, and matched the trailing `135` as a separate
  value — which is how a literature value appeared to be present in a sentence that never
  mentioned it.
- **Only numbers above 1000 were required to be traceable.** "20 protein-coding genes",
  "25 tRNAs", "6 introns", "41 modules" were never checked.
- **A number merely had to exist somewhere in the fact dictionary.** Nothing tied a number to
  the claim it was making.

All fixed. The threshold is gone: every number must now be a computed fact, a declared
parameter, a value attributed to a cited paper, or an explicitly reviewed value with recorded
provenance.

## 10. New and rebuilt checks

| Script | What it enforces |
|---|---|
| `scripts/50_audit_citations.py` | Every DOI resolves at CrossRef and matches on title, full author list, year, volume and first page. First-author mismatch is called out as the signature of an entry attached to the wrong paper. Handles `and others`, corporate authors and LaTeX accents; tolerates CrossRef records truncated to one author rather than failing on them. Also checks every `\cite` key resolves and every entry is used. |
| `scripts/51_audit_identifiers.py` | NCBI accessions against local FASTA headers, then the nuccore or datasets API, with an organism check; Rfam IDs against the `DESC` line of the downloaded covariance models; EC numbers against the Swiss-Prot flat file; ModelSEED IDs against `refs/modelseed/`. Undoes LaTeX escaping first, without which `GCA\_056149245.1` is missed. |
| `scripts/32_audit_numbers.py` | 312 facts computed live from result files. Adds (a) **claim-attachment assertions** — 25 headline numbers must appear in the sentence making the claim, not merely in the document; (b) an **EXTERNAL** category requiring each literature value to sit in a sentence citing its source; (c) a **REVIEWED** allowlist recording provenance for hand-verified values; (d) a **cross-document consistency** check comparing the 57 sentences shared between manuscripts with their numbers removed. |

Run all three before any submission build. Each exits non-zero on failure.

## 11. Verified and unchanged

Checked and found correct, listed because a reader of this report is entitled to know the
denominator:

- All six EC numbers, all five Rfam IDs, all seven NCBI accessions — including that each
  assembly accession belongs to the strain the text names for it.
- The mito:nuclear expression ratios (median 1.51 exudophore; 0.53–0.74 elsewhere; 1.02, 0.94
  and 1.14 for the individual libraries called out) — all reproduce as per-library medians.
- The comparative correlations `r = +0.89` (intron content vs length) and `+0.96` (ORF count vs
  length), the 1.7-fold length range, the 7.8 kb within-species spread, 26.6% GC.
- 42.0–96.6% rRNA of assigned counts; 82 → 281 genes with a significant tissue effect on
  filtering; 92,430,743 raw reads; 2,135-fold assignment gain; 53 features / 47 quantifiable.
- Causal language: the manuscripts are already appropriately hedged, and P5 explicitly declines
  to claim retrograde signalling from expression balance.
- Replicate transparency: P2, P3 and P5 each state that the exudophore rests on two libraries
  after filtering.
- Cross-document consistency: 57 sentences shared between manuscripts, no quantity stated two
  ways.

## 12. Still outstanding — needs the authors, not the audit

Unchanged by this work and still blocking submission: exudophore morphology and the sampling
protocol; RNA extraction method; culture age at sampling; ORCIDs; affiliations; CRediT roles;
funding; competing interests. The given-name/surname splits for Leiva, Dagar and Rizwan remain
unconfirmed (those for Sanchez and Moulton were resolved — see §3).
