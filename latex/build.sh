#!/usr/bin/env bash
# Build any manuscript target. Usage: build.sh [v0|p1|p2|all]
#
# refs.bib and the .bst files are copied into the target directory rather than referenced by
# relative path: bibtex resolves \bibliography{} relative to its own working directory, so
# \bibliography{../refs} silently produces an empty bibliography.
#
# No associative arrays: macOS ships bash 3.2, which does not support them.
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
TEC="$ROOT/../envs/root/envs/env-tex/bin/tectonic"

dir_for()  { case "$1" in v0) echo V0_internal;; p1) echo P1_resource;; p2) echo P2_biology;; p3) echo P3_systems;; p4) echo P4_mitochondrial;; p5) echo P5_mitocomparative;; esac; }
name_for() { case "$1" in
  v0) echo Myco_RNAseq_internal_draft;;
  p1) echo Myco_RNAseq_resource_draft;;
  p2) echo Myco_exudophore_biology_draft;;
  p3) echo Myco_systems_metabolism_draft;;
  p4) echo Myco_mitogenome_draft;;
  p5) echo Myco_mitocomparative_draft;; esac; }

build_one() {
  key="$1"; d="$ROOT/$(dir_for "$key")"; out="$(name_for "$key")"
  if [ ! -f "$d/main.tex" ]; then echo "[skip] $key (no main.tex yet)"; return 0; fi
  echo "[build] $key"
  cp "$ROOT/refs.bib" "$d/"
  cp "$ROOT/template/sn-jnl.cls" "$d/"
  cp "$ROOT/template/bst/"*.bst "$d/"
  ( cd "$d" && "$TEC" -X compile main.tex --keep-intermediates --keep-logs >/dev/null 2>&1 || true )
  ( cd "$d" && "$TEC" -X compile main.tex --keep-intermediates --keep-logs >/dev/null 2>&1 )
  und=$(grep -c "undefined" "$d/main.log" 2>/dev/null | tr -d ' ' || true)
  [ "${und:-0}" != "0" ] && echo "  note: ${und} undefined reference/citation warning(s)"
  cites=$(grep -c "bibitem" "$d/main.bbl" 2>/dev/null | tr -d ' ' || true)
  echo "  bibliography entries: ${cites:-0}"
  mkdir -p "$ROOT/../manuscript"
  cp "$d/main.pdf" "$ROOT/../manuscript/${out}.pdf"
  echo "  -> manuscript/${out}.pdf ($(du -h "$d/main.pdf" | cut -f1 | tr -d ' '))"
}

# The three audits gate the build. `audit` runs them alone; `all` runs them first and stops
# if any fails, so a submission PDF cannot be produced from sources with an unverified
# citation, identifier or number. Use `all-noaudit` to build while still editing.
run_audits() {
  local rc=0
  for s in 50_audit_citations 51_audit_identifiers 32_audit_numbers; do
    echo "[audit] $s"
    python3 "$ROOT/../scripts/${s}.py" > "/tmp/${s}.out" 2>&1 || { rc=1; }
    tail -1 "/tmp/${s}.out" | sed 's/^/  /'
    [ "$rc" = 0 ] || { echo "  FAILED - see /tmp/${s}.out"; return 1; }
  done
  return 0
}

case "${1:-all}" in
  audit)      run_audits ;;
  all)        run_audits || { echo "audit failed; not building"; exit 1; }
              for k in v0 p1 p2 p3 p4 p5; do build_one "$k"; done ;;
  all-noaudit) for k in v0 p1 p2 p3 p4 p5; do build_one "$k"; done ;;
  *)          build_one "${1}" ;;
esac
