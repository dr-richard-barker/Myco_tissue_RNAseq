#!/usr/bin/env bash
# Phase 3 -- full pipeline over all 16 samples, adapted from NASA GeneLab GL-DPPD-7101-G.
#
# Usage: 10_pipeline.sh [reference_label]      (default PC9.15; the other is BOM_ss5)
#
# Trimming is reference-independent, so trimmed FASTQs in qc/fastp/ are shared and reused
# across references; alignment onwards is written per reference into bam_<label>/.
#
# Documented deviations from the GeneLab spec (rationale in NOTES.md):
#   * HISAT2 replaces STAR -- STAR is non-functional on this machine (see 02_testmap.sh).
#   * featureCounts replaces RSEM  -- RSEM models full-length transcript coverage and is
#     invalid for 3'-tag pileups.
#   * fastp replaces TrimGalore!   -- the vendor already adapter-trimmed; the real problems
#     are polyG (NovaSeq two-colour) and polyA read-through, which fastp handles directly.
#   * counting annotation is <label> + 3' extension + barrnap/coverage-derived rRNA features.
#   * -s 1 (forward-stranded), established empirically: -s1 49.6% vs -s2 6.3% assigned.
#
# Envs go on PATH rather than being invoked via `micromamba run` per command: the latter
# takes a package-cache lock on every call, which serialised the whole run.
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
export PATH="$ROOT/envs/root/envs/env-hisat/bin:$PATH"
UMI="$ROOT/envs/root/envs/env-qc/bin/umi_tools"

LABEL="${1:-PC9.15}"
FASTQ_DIR="$ROOT/../GPNJ7M_fastq"
GTF="$ROOT/refs/$LABEL/${LABEL}_final.gtf"
IDX="$ROOT/refs/$LABEL/ht2"
BAMDIR="$ROOT/bam_$LABEL"
THREADS=14

[[ -s "$GTF" ]] || { echo "missing $GTF -- run 08_build_annotation.py first"; exit 1; }
mkdir -p "$ROOT"/{qc/fastp,counts,logs} "$BAMDIR"
echo "=== reference: $LABEL ==="

for fq in "$FASTQ_DIR"/*.fastq.gz; do
  s=$(basename "$fq" .fastq.gz)
  trimmed="$ROOT/qc/fastp/$s.trim.fq.gz"
  part="$ROOT/qc/fastp/$s.part.fq.gz"   # must still end in .gz: fastp picks compression
  bam="$BAMDIR/$s.bam"                   # from the filename extension, not a flag
  prim="$BAMDIR/$s.primary.bam"
  dedup="$BAMDIR/$s.dedup.bam"

  if [[ ! -s "$trimmed" ]]; then
    echo "[trim ] $s"
    fastp -i "$fq" -o "$part" \
      --trim_poly_g --trim_poly_x --poly_x_min_len 8 \
      --cut_tail --cut_tail_mean_quality 15 --length_required 30 \
      --json "$ROOT/qc/fastp/$s.json" --html "$ROOT/qc/fastp/$s.html" \
      --thread 8 2> "$ROOT/logs/fastp_$s.log"
    gzip -t "$part"   # fail loudly rather than feeding a plain file to the aligner
    mv "$part" "$trimmed"
  fi

  if [[ ! -s "$bam" ]]; then
    echo "[align] $s"
    hisat2 -p $THREADS -x "$IDX" -U "$trimmed" \
      --max-intronlen 3000 --rna-strandness F \
      --new-summary --summary-file "$ROOT/qc/${LABEL}__$s.hisat2.txt" 2>> "$ROOT/logs/hisat2_$LABEL.log" \
      | samtools sort -@ 4 -o "$bam.part" - 2>> "$ROOT/logs/hisat2_$LABEL.log"
    mv "$bam.part" "$bam"
    samtools index "$bam"
    n=$(samtools view -c -F 4 "$bam")
    [[ "$n" -gt 0 ]] || { echo "[FATAL] $s aligned 0 reads"; exit 1; }
  fi

  # Keep only primary alignments before dedup. HISAT2 emits ~8.9M secondary records per
  # 2.5M primary ones here, almost entirely rRNA multi-mapping across repeated rDNA copies.
  # featureCounts discards multi-mappers regardless, and feeding them to umi_tools made
  # dedup ~5x slower (a projected 3.7 h across the run).
  if [[ ! -s "$prim" ]]; then
    samtools view -b -F 0x104 -@ 4 "$bam" -o "$prim.part"
    mv "$prim.part" "$prim"
    samtools index "$prim"
  fi

  # UMI deduplication. The 14 nt UMI is appended to the read name after '_' by the vendor.
  if [[ ! -s "$dedup" ]]; then
    echo "[dedup] $s"
    if "$UMI" dedup --stdin="$prim" --stdout="$dedup.part" \
         --umi-separator=_ --method=directional \
         > "$ROOT/logs/umitools_${LABEL}_$s.log" 2>&1; then
      mv "$dedup.part" "$dedup"
      samtools index "$dedup"
    else
      echo "[WARN ] umi_tools failed for $s -- see logs/umitools_${LABEL}_$s.log"
      rm -f "$dedup.part"
    fi
  fi
done

for mode in dedup raw; do
  if [[ "$mode" == dedup ]]; then
    bams=("$BAMDIR"/*.dedup.bam)
  else
    bams=("$BAMDIR"/*.primary.bam)
  fi
  [[ -e "${bams[0]}" ]] || continue
  echo "[count] $mode (${#bams[@]} bams)"
  featureCounts -T $THREADS -t exon -g gene_id -s 1 \
    -a "$GTF" -o "$ROOT/counts/counts_${LABEL}_$mode.txt" "${bams[@]}" \
    > "$ROOT/logs/featurecounts_${LABEL}_$mode.log" 2>&1
done

echo "=== pipeline complete ($LABEL) ==="
