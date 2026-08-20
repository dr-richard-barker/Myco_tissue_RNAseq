#!/usr/bin/env bash
# P5 -- fetch a comparative set of complete Pleurotus mitogenomes.
#
# BOM_ss14 is included deliberately as a control: it is the sibling nucleus of the same
# dikaryon as BOM_ss5 (71,947 vs 71,949 bp), so it should be near-identical. Material
# divergence between them would mean the comparative pipeline is broken rather than that the
# nuclei differ.
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
OUT="$ROOT/refs/mito_comparative"; mkdir -p "$OUT"
E="https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

# Verified Pleurotus mitogenome accessions only. Two accessions in an earlier draft of this
# list (NC_037259.1, NC_042144.1) were guessed rather than looked up and turned out to be a
# Physcomitrella chromosome and a Cymbopogon chloroplast; every accession here has been
# checked against its FASTA header.
ACC="CM148778.1 PX724300.1 PX724301.1 PX724302.1 OR030114.1 NC_062374.1 NC_061177.1 \
NC_038091.1 NC_036999.1 NC_036998.1"

for a in $ACC; do
  f="$OUT/$a.fa"
  [[ -s "$f" ]] && { echo "[skip] $a"; continue; }
  curl -s -G "$E/efetch.fcgi" --data-urlencode "db=nuccore" --data-urlencode "id=$a" \
    -d rettype=fasta -d retmode=text -o "$f.part" --max-time 120 || true
  if [[ -s "$f.part" ]] && head -1 "$f.part" | grep -q "^>"; then
    mv "$f.part" "$f"; echo "[get ] $a"
  else
    rm -f "$f.part"; echo "[MISS] $a"
  fi
  sleep 0.4
done
# our own two, from the local assemblies
cp "$ROOT/refs/mito/CM148777.1.fa" "$OUT/CM148777.1.fa" 2>/dev/null || true
python3 - <<'PY'
import pathlib
out = pathlib.Path("/Users/drb_laptop/Documents/Tissue_specific_myeclium/analysis/refs/mito_comparative")
seqs={}; n=None
for l in open("/Users/drb_laptop/Documents/Tissue_specific_myeclium/analysis/refs/PC9.15/PC9.15_genomic.fna"):
    if l[0]==">": n=l[1:].split()[0]; seqs[n]=[]
    elif n: seqs[n].append(l.strip())
s="".join(seqs.get("CM057219.1",[]))
if s:
    (out/"CM057219.1.fa").write_text(">CM057219.1 Pleurotus ostreatus PC9.15 mitochondrion\n" +
        "\n".join(s[i:i+70] for i in range(0,len(s),70)) + "\n")
PY
echo "=== set ==="
# guard: refuse anything that is not a Pleurotus mitogenome of plausible size
for f in "$OUT"/*.fa; do
  n=$(grep -v '^>' "$f" | tr -d '\n' | wc -c | tr -d ' ')
  if ! head -1 "$f" | grep -qi "pleurotus" || [ "$n" -gt 200000 ]; then
    echo "[DROP] $(basename "$f" .fa) -- not a Pleurotus mitogenome ($n bp)"; rm -f "$f"
  fi
done
for f in "$OUT"/*.fa; do
  n=$(grep -v '^>' "$f" | tr -d '\n' | wc -c | tr -d ' ')
  t=$(head -1 "$f" | cut -c2-62)
  printf "  %-14s %9s bp  %s\n" "$(basename "$f" .fa)" "$n" "$t"
done
