#!/usr/bin/env bash
# Phase 0 -- build conda environments.
# Tool versions pinned to NASA GeneLab GL-DPPD-7101-G where that pipeline specifies them.
# Most bioconda RNA-seq tools have no osx-arm64 build, so the align/QC envs are built as
# osx-64 and run under Rosetta 2 (verified available on this machine).
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
export MAMBA_ROOT_PREFIX="$ROOT/envs/root"
MM="$ROOT/envs/bin/micromamba"
CH="-c conda-forge -c bioconda"

# Compiled aligner/counting tools. barrnap lives here too (perl + nhmmer, no python pins).
"$MM" create -y -n env-align --platform osx-64 $CH \
  star=2.7.11b samtools=1.21 subread=2.1.1 bedtools seqkit \
  ucsc-gtftogenepred ucsc-genepredtobed barrnap

# Python/Java QC + trimming tools, kept separate so their python pins cannot
# conflict with the solver for env-align.
"$MM" create -y -n env-qc --platform osx-64 $CH \
  fastqc=0.12.1 multiqc rseqc=5.0.4 qualimap cutadapt=4.2 umi_tools

# R stack, native arm64.
"$MM" create -y -n env-r $CH \
  r-base=4.4 bioconductor-deseq2 bioconductor-edger bioconductor-tximport \
  r-tidyverse r-pheatmap r-ggrepel

echo "=== environments built ==="
"$MM" env list
