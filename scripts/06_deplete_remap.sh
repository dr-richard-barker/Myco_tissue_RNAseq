#!/usr/bin/env bash
# Phase 1/2a -- measure (not estimate) the unique mapping rate after rRNA depletion.
#
# The Phase 1 gate of ">60% uniquely mapped" was written assuming an ordinary RNA-seq
# library. These libraries are ~53% rRNA, and rDNA is a tandem repeat, so rRNA reads pile
# into the multi-mapping bin and crush the raw unique rate. The meaningful question is what
# fraction of the NON-rRNA reads map uniquely, so deplete first and re-map.
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
export MAMBA_ROOT_PREFIX="$ROOT/envs/root"
MM="$ROOT/envs/bin/micromamba"
LABEL="${1:-BOM_ss5}"

SUB="$ROOT/qc/subsample"
DEP="$ROOT/qc/depleted"
OUT="$ROOT/qc/testmap_depleted"
mkdir -p "$DEP" "$OUT"

rrna_idx="$ROOT/refs/rRNA/${LABEL}_rRNA_ht2"
gen_idx="$ROOT/refs/$LABEL/ht2"

for fq in "$SUB"/*.fq; do
  s=$(basename "$fq" .fq)
  if [[ ! -s "$DEP/$s.fq" ]]; then
    "$MM" run -n env-hisat hisat2 -p 14 -x "$rrna_idx" -U "$fq" \
      --un "$DEP/$s.fq.part" -S /dev/null \
      --summary-file "$ROOT/qc/rrna/${s}.summary.txt" --new-summary \
      >> "$ROOT/logs/ht2_deplete.log" 2>&1
    mv "$DEP/$s.fq.part" "$DEP/$s.fq"
  fi
  sm="$OUT/${LABEL}__${s}.summary.txt"
  [[ -s "$sm" ]] && continue
  "$MM" run -n env-hisat hisat2 -p 14 -x "$gen_idx" -U "$DEP/$s.fq" \
    --max-intronlen 3000 --new-summary --summary-file "$sm" -S /dev/null \
    >> "$ROOT/logs/ht2_deplete.log" 2>&1
done

echo "=== depleted remap done for $LABEL ==="
