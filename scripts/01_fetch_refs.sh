#!/usr/bin/env bash
# Phase 1 -- fetch the four candidate Pleurotus reference genomes.
# The species call from the reads is genus-level only (28S LSU is 100% identical across
# Pleurotus), so the reference is chosen empirically on mapping rate, not assumed.
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT/refs"

# accession                 label      ftp_dir_path
CANDIDATES=(
  "GCF_014466165.1 PC9      GCF/014/466/165/GCF_014466165.1_ASM1446616v1"
  "GCA_029852705.2 PC9.15   GCA/029/852/705/GCA_029852705.2_AS_PC9.15"
  "GCA_056149245.1 BOM_ss5  GCA/056/149/245/GCA_056149245.1_AS_BOM_ss5"
  "GCA_056149315.1 BOM_ss14 GCA/056/149/315/GCA_056149315.1_AS_BOM_ss14"
)

BASE="https://ftp.ncbi.nlm.nih.gov/genomes/all"

for entry in "${CANDIDATES[@]}"; do
  read -r acc label path <<< "$entry"
  stem="$(basename "$path")"
  mkdir -p "$label"
  for ext in genomic.fna.gz genomic.gtf.gz protein.faa.gz; do
    out="$label/${label}_${ext}"
    if [[ -s "$out" ]]; then
      echo "[skip] $out"
      continue
    fi
    url="$BASE/$path/${stem}_${ext}"
    echo "[get ] $label $ext"
    if ! curl -fsSL --retry 3 -o "$out.part" "$url"; then
      echo "[MISS] $label $ext not available at $url"
      rm -f "$out.part"
      continue
    fi
    mv "$out.part" "$out"
  done
done

echo
echo "=== fetched ==="
find . -name '*.gz' -type f | sort | while read -r f; do
  printf "%-52s %s\n" "$f" "$(du -h "$f" | cut -f1)"
done
