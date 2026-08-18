#!/usr/bin/env bash
# Prepare derived tables that the figure script needs but that no earlier step wrote to disk.
# Kept separate from 30_figures.R so the (slow) BAM traversals run once and the figures can be
# re-rendered instantly from committed CSVs.
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
export PATH="$ROOT/envs/root/envs/env-hisat/bin:$PATH"
OUT="$ROOT/results/figure_data"; mkdir -p "$OUT"

# --- 1. per-base coverage over the BOM_ss5 mitochondrion (Fig 2c) -------------------------
# The mitochondrial rDNA block is unannotated; this trace is the evidence for it.
if [[ ! -s "$OUT/mito_coverage.csv" ]]; then
  echo "[cov] mitochondrial coverage"
  # pool three representative libraries spanning the yield range
  samtools merge -f -@ 4 /tmp/mito_pool.bam \
    "$ROOT"/bam_BOM_ss5/GPNJ7M_16_sample_2H.primary.bam \
    "$ROOT"/bam_BOM_ss5/GPNJ7M_1_sample_1A.primary.bam \
    "$ROOT"/bam_BOM_ss5/GPNJ7M_11_sample_2C.primary.bam 2>/dev/null
  samtools index /tmp/mito_pool.bam
  { echo "pos,depth"
    samtools depth -a -r CM148777.1 /tmp/mito_pool.bam | awk '{print $2","$3}'
  } > "$OUT/mito_coverage.csv"
  rm -f /tmp/mito_pool.bam /tmp/mito_pool.bam.bai
fi

# --- 2. annotated features on the mitochondrion, for the track beneath (Fig 2c) -----------
if [[ ! -s "$OUT/mito_features.csv" ]]; then
  { echo "start,end,strand,type,name"
    awk -F'\t' '$1=="CM148777.1" && $3=="gene"{
      bt="other"; if (match($9,/gene_biotype "[^"]+"/)) bt=substr($9,RSTART+15,RLENGTH-16)
      nm="";      if (match($9,/gene_id "[^"]+"/))      nm=substr($9,RSTART+9,RLENGTH-10)
      print $4","$5","$7","bt","nm }' "$ROOT/refs/BOM_ss5/BOM_ss5_genomic.gtf"
  } > "$OUT/mito_features.csv"
fi

# --- 3. UMI duplication per library (Supp Fig S6) -----------------------------------------
if [[ ! -s "$OUT/dedup_rates.csv" ]]; then
  echo "[dup] duplication rates"
  { echo "sample,primary,dedup"
    for f in "$ROOT"/bam_BOM_ss5/*.dedup.bam; do
      s=$(basename "$f" .dedup.bam)
      printf "%s,%s,%s\n" "$s" "$(samtools view -c "$ROOT/bam_BOM_ss5/$s.primary.bam")" "$(samtools view -c "$f")"
    done
  } > "$OUT/dedup_rates.csv"
fi

# --- 4. intergenic peaks that revealed the rDNA blocks (Supp Fig S7) -----------------------
if [[ ! -s "$OUT/intergenic_peaks.csv" ]]; then
  echo "[peak] intergenic clustering"
  awk -F'\t' '$3=="gene"{print $1"\t"($4-1)"\t"$5}' "$ROOT/refs/BOM_ss5/BOM_ss5_genomic.gtf" \
    | sort -k1,1 -k2,2n > /tmp/genes.bed
  samtools merge -f -@ 4 /tmp/peak_pool.bam "$ROOT"/bam_BOM_ss5/*.primary.bam 2>/dev/null
  TOT=$(samtools view -c -F 4 /tmp/peak_pool.bam)
  "$ROOT/envs/root/envs/env-align/bin/bedtools" intersect -v -abam /tmp/peak_pool.bam -b /tmp/genes.bed 2>/dev/null \
    | samtools view -F 4 - 2>/dev/null | awk '{print $3"\t"$4"\t"$4+length($10)}' | sort -k1,1 -k2,2n > /tmp/ig.bed
  { echo "contig,start,end,reads,pct_of_aligned"
    "$ROOT/envs/root/envs/env-align/bin/bedtools" merge -d 300 -c 1 -o count -i /tmp/ig.bed \
      | awk -v T="$TOT" '$4>=500{printf "%s,%s,%s,%s,%.4f\n",$1,$2,$3,$4,100*$4/T}' | sort -t, -k4 -rn
  } > "$OUT/intergenic_peaks.csv"
  rm -f /tmp/peak_pool.bam /tmp/ig.bed /tmp/genes.bed
fi

echo "=== figure data ready ==="
wc -l "$OUT"/*.csv
