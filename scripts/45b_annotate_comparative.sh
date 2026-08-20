#!/usr/bin/env bash
# Annotate every comparative mitogenome with OUR pipeline rather than trusting deposited
# features. Deposited annotation quality varies enormously -- that is P4's whole finding -- so
# comparing our annotation against theirs would confound biology with annotation effort.
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
export MAMBA_ROOT_PREFIX="$ROOT/envs/root"
MM="$ROOT/envs/bin/micromamba"
IN="$ROOT/refs/mito_comparative"; OUT="$ROOT/results/mito_comparative"; mkdir -p "$OUT"

echo "acc,length,gc_pct,tRNA,orfs_ge100aa,orfs_with_hit" > "$OUT/inventory.csv"
for f in "$IN"/*.fa; do
  acc=$(basename "$f" .fa)
  len=$(grep -v '^>' "$f" | tr -d '\n' | wc -c | tr -d ' ')
  gc=$(python3 -c "
s=''.join(l.strip() for l in open('$f') if not l.startswith('>')).upper()
print(f'{100*(s.count(\"G\")+s.count(\"C\"))/len(s):.1f}')")
  # tRNAs
  if [[ ! -s "$OUT/$acc.trna" ]]; then
    "$MM" run -n env-mito tRNAscan-SE -O -q -o "$OUT/$acc.trna" "$f" > /dev/null 2>&1 || true
  fi
  nt=$(grep -c "$acc" "$OUT/$acc.trna" 2>/dev/null | tr -d ' ' || echo 0)
  # ORFs under genetic code 4
  if [[ ! -s "$OUT/$acc.orfs.fa" ]]; then
    "$MM" run -n env-mito getorf -sequence "$f" -table 4 -minsize 300 -find 1 \
      -outseq "$OUT/$acc.orfs.fa" 2>/dev/null || true
  fi
  no=$(grep -c '^>' "$OUT/$acc.orfs.fa" 2>/dev/null | tr -d ' ' || echo 0)
  # identify them
  if [[ ! -s "$OUT/$acc.hits.tsv" && "$no" -gt 0 ]]; then
    "$MM" run -n env-annot diamond blastp -q "$OUT/$acc.orfs.fa" -d "$ROOT/refs/uniprot/sprot" \
      -o "$OUT/$acc.hits.tsv" --threads 8 --max-target-seqs 1 --evalue 1e-5 --quiet \
      --outfmt 6 qseqid sseqid pident length evalue bitscore stitle > /dev/null 2>&1 || true
  fi
  nh=$(cut -f1 "$OUT/$acc.hits.tsv" 2>/dev/null | sort -u | wc -l | tr -d ' ' || echo 0)
  echo "$acc,$len,$gc,$nt,$no,$nh" >> "$OUT/inventory.csv"
  printf "  %-14s %7s bp  GC %s%%  tRNA %-3s ORFs %-3s identified %s\n" "$acc" "$len" "$gc" "$nt" "$no" "$nh"
done
