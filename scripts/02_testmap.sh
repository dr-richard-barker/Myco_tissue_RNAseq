#!/usr/bin/env bash
# Phase 1 -- decide the reference on mapping rate rather than assumption.
#
# ALIGNER DEVIATION FROM GL-DPPD-7101-G: this uses HISAT2, not STAR.
# STAR is non-functional on this machine (macOS Darwin 25.5 / Apple Silicon). Both the
# bioconda osx-64 build (under Rosetta) and the native osx-arm64 build load the genome
# correctly but then report "Number of input reads | 0" with "nextChar=-1" for every input,
# including a synthetic 200-read/200 kb genome control, at --runThreadN 1 and 2, with and
# without --readFilesCommand, and with the sandbox disabled. seqkit in the same env reads
# the identical files without trouble, so it is STAR itself, not the environment.
# HISAT2 2.2.3 passes the same synthetic control at 100% alignment. It is splice-aware and
# soft-clips by default, which Phase 2 needs to call polyadenylation sites.
#
# Genome-only indices so this measures genome fit alone, not annotation quality (which
# differs sharply between candidates: PC9.15 has real UTRs, the other three do not).
# 250k reads/sample; the same reads go to every reference so head-of-file tile bias cancels.
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
export MAMBA_ROOT_PREFIX="$ROOT/envs/root"
MM="$ROOT/envs/bin/micromamba"

FASTQ_DIR="$ROOT/../GPNJ7M_fastq"
SUB="$ROOT/qc/subsample"
OUT="$ROOT/qc/testmap"
NREADS=250000
THREADS=14

mkdir -p "$SUB" "$OUT"

for fq in "$FASTQ_DIR"/*.fastq.gz; do
  s=$(basename "$fq" .fastq.gz)
  [[ -s "$SUB/$s.fq" ]] && continue
  # head closes the pipe early, so gzcat takes SIGPIPE; expected, not an error.
  { gzcat "$fq" || true; } | head -n $((NREADS * 4)) > "$SUB/$s.fq.part"
  mv "$SUB/$s.fq.part" "$SUB/$s.fq"
done
echo "subsampled $(ls "$SUB"/*.fq | wc -l | tr -d ' ') samples to ${NREADS} reads"

for label in PC9 PC9.15 BOM_ss5 BOM_ss14; do
  idx="$ROOT/refs/$label/ht2"
  fna="$ROOT/refs/$label/${label}_genomic.fna"
  [[ -s "$fna" ]] || gzcat "$fna.gz" > "$fna"
  if [[ ! -s "$idx.1.ht2" ]]; then
    echo "[index] $label"
    "$MM" run -n env-hisat hisat2-build -p $THREADS "$fna" "$idx" \
      > "$ROOT/logs/ht2_index_$label.log" 2>&1
  fi
  for fq in "$SUB"/*.fq; do
    s=$(basename "$fq" .fq)
    sm="$OUT/${label}__${s}.summary.txt"
    [[ -s "$sm" ]] && continue
    "$MM" run -n env-hisat hisat2 -p $THREADS -x "$idx" -U "$fq" \
      --max-intronlen 3000 --new-summary --summary-file "$sm" -S /dev/null \
      >> "$ROOT/logs/ht2_testmap.log" 2>&1
  done
  echo "[mapped] $label"
done

echo "=== done; parse with 03_testmap_report.py ==="
