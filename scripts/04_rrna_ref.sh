#!/usr/bin/env bash
# Phase 2a -- locate rRNA loci and measure the true rRNA fraction per sample.
#
# Needed because the test-map shows ~60% multi-mapping reads: rDNA is a tandem repeat, so
# rRNA reads hit many copies and land in the multi bin. The probe-based floor computed from
# the raw FASTQs (12-46%) counted only four diagnostic 25-mers and was always a lower bound.
# This replaces it with a direct measurement, and produces the rRNA locus list that Phase 4
# needs for the GL-DPPD-7101-G rRNA-removed DGE track.
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
export MAMBA_ROOT_PREFIX="$ROOT/envs/root"
MM="$ROOT/envs/bin/micromamba"
LABEL="${1:?usage: 04_rrna_ref.sh <reference-label>}"

REF="$ROOT/refs/$LABEL"
FNA="$REF/${LABEL}_genomic.fna"
OUT="$ROOT/refs/rRNA"
SUB="$ROOT/qc/subsample"
mkdir -p "$OUT"

[[ -s "$FNA" ]] || gzcat "$FNA.gz" > "$FNA"

# barrnap euk finds 18S/5.8S/28S/5S; mito finds the organellar rRNAs. Both matter: the
# mitochondrion is a separate high-copy compartment with its own rRNA.
for kingdom in euk mito; do
  gff="$OUT/${LABEL}_${kingdom}.gff"
  [[ -s "$gff" ]] && continue
  "$MM" run -n env-align barrnap --kingdom "$kingdom" --threads 8 "$FNA" \
    > "$gff" 2> "$ROOT/logs/barrnap_${LABEL}_${kingdom}.log"
done

cat "$OUT/${LABEL}_euk.gff" "$OUT/${LABEL}_mito.gff" | grep -v '^#' | sort -k1,1 -k4,4n \
  > "$OUT/${LABEL}_rRNA.gff" || true

awk -F'\t' '!/^#/ && NF>=9 {
  name="rRNA"; if (match($9, /Name=[^;]+/)) name=substr($9, RSTART+5, RLENGTH-5)
  print $1"\t"($4-1)"\t"$5"\t"name"\t0\t"$7
}' "$OUT/${LABEL}_rRNA.gff" | sort -k1,1 -k2,2n > "$OUT/${LABEL}_rRNA.bed"

echo "=== rRNA features found in $LABEL ==="
cut -f4 "$OUT/${LABEL}_rRNA.bed" | sort | uniq -c | sort -rn

"$MM" run -n env-align bedtools getfasta -s -fi "$FNA" \
  -bed "$OUT/${LABEL}_rRNA.bed" -fo "$OUT/${LABEL}_rRNA.fa" -name

idx="$OUT/${LABEL}_rRNA_ht2"
[[ -s "$idx.1.ht2" ]] || "$MM" run -n env-hisat hisat2-build -p 8 \
  "$OUT/${LABEL}_rRNA.fa" "$idx" > "$ROOT/logs/ht2_index_rRNA_$LABEL.log" 2>&1

mkdir -p "$ROOT/qc/rrna"
for fq in "$SUB"/*.fq; do
  s=$(basename "$fq" .fq)
  sm="$ROOT/qc/rrna/${s}.summary.txt"
  [[ -s "$sm" ]] && continue
  "$MM" run -n env-hisat hisat2 -p 14 -x "$idx" -U "$fq" \
    --new-summary --summary-file "$sm" -S /dev/null \
    >> "$ROOT/logs/ht2_rrna.log" 2>&1
done

echo "=== done; parse with 05_rrna_report.py ==="
