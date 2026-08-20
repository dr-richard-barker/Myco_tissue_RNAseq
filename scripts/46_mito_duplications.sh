#!/usr/bin/env bash
# P5 -- internal repeats/duplications, and pairwise synteny against the comparative set.
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
export MAMBA_ROOT_PREFIX="$ROOT/envs/root"
MM="$ROOT/envs/bin/micromamba"
REF="$ROOT/refs/mito/CM148777.1.fa"
OUT="$ROOT/results/mito_dup"; mkdir -p "$OUT"

# --- self-alignment: internal repeats. --maxmatch keeps non-unique matches, which is the
# whole point here; the trivial full-length diagonal is filtered out downstream.
"$MM" run -n env-mito nucmer --maxmatch --nosimplify -p "$OUT/self" "$REF" "$REF" 2>/dev/null
"$MM" run -n env-mito show-coords -rclTH "$OUT/self.delta" > "$OUT/self.coords" 2>/dev/null
echo "  self-alignments (pre-filter): $(wc -l < "$OUT/self.coords" | tr -d ' ')"

# --- pairwise against every comparative genome
for f in "$ROOT"/refs/mito_comparative/*.fa; do
  acc=$(basename "$f" .fa); [ "$acc" = "CM148777.1" ] && continue
  "$MM" run -n env-mito nucmer -p "$OUT/vs_$acc" "$REF" "$f" 2>/dev/null
  "$MM" run -n env-mito show-coords -rclTH "$OUT/vs_$acc.delta" > "$OUT/vs_$acc.coords" 2>/dev/null
done
echo "  pairwise alignments done: $(ls "$OUT"/vs_*.coords 2>/dev/null | wc -l | tr -d ' ')"
