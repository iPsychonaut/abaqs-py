#!/usr/bin/env bash
#
# Differential parity run: execute the Java ABAQS and the Python port on the
# EXACT same inputs, capture both sets of metrics into one combined report,
# and run the two-tier parity check (strict display + numeric tolerance).
#
# Prereq: run setup_env.sh once, then `conda activate abaqs`.
#
# Put the downloaded JGI files in $RAW_DIR (or point RAW_DIR at them), e.g.:
#   RAW_DIR=/mnt/c/Users/theda/Downloads bash run_parity.sh
#
set -euo pipefail

# ---------- Configuration (override via env) ----------
ABAQS_JAVA_DIR="${ABAQS_JAVA_DIR:-/mnt/c/Users/theda/Desktop/abaqs-java}"
JAR="${JAR:-$ABAQS_JAVA_DIR/target/abaqs-jar-with-dependencies.jar}"
PY_ABAQS="${PY_ABAQS:-abaqs}"                 # console script from `pip install -e`
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
COMPARE="${COMPARE:-$SCRIPT_DIR/compare_abaqs.py}"

RAW_DIR="${RAW_DIR:-$HOME/abaqs-parity/raw}"   # where the downloaded JGI files live
OUT_DIR="${OUT_DIR:-$HOME/abaqs-parity/out}"   # where results are written
GENOME="${GENOME:-Altalt1}"
TOL="${TOL:-0.0001}"                           # Tier-2 tolerance = 1 unit of last displayed digit

# Fixed BUSCO string handed IDENTICALLY to both tools (parity, not biology).
BUSCO_STRING="${BUSCO_STRING:-C:99.3%[S:98.9%,D:0.4%],F:0.3%,M:0.4%,n:758}"

mkdir -p "$OUT_DIR"
[[ -f "$JAR" ]] || { echo "ERROR: jar not found: $JAR  (run setup_env.sh)" >&2; exit 2; }

# ---------- Resolve inputs ----------
# Find the first file matching a (case-insensitive) glob anywhere under RAW_DIR.
# Recursive so it works whether the Globus download is flat or keeps its nested
# Altalt1/Mycocosm/Annotation/... directory layout.
find_one() { find "$RAW_DIR" -type f -iname "$1" 2>/dev/null | sort | head -n1; }

# GFF3: prefer the clean, FILTERED GeneCatalog .gff3.gz (both tools read .gz
# directly). Fall back to an extracted/tarballed all-models GFF3 if that's all
# that's present. NB: never the *.gff.gz — that's the old GFF2 format.
GFF3="${GFF3:-}"
if [[ -z "$GFF3" ]]; then
  GFF3="$(find_one "${GENOME}_GeneCatalog_*.gff3.gz")"
  [[ -z "$GFF3" ]] && GFF3="$(find_one "${GENOME}_GeneCatalog_*.gff3")"
  if [[ -z "$GFF3" ]]; then
    tgz="$(find_one "${GENOME}_all_genes_*.gff3.tgz")"
    if [[ -n "$tgz" ]]; then
      echo "==> Extracting GFF3 from tarball…"
      tar -xzf "$tgz" -C "$RAW_DIR"
      GFF3="$(find "$RAW_DIR" -type f -iname '*.gff3' | sort | head -n1)"
    fi
  fi
fi
[[ -n "$GFF3" && -f "$GFF3" ]] || { echo "ERROR: GFF3 not found under $RAW_DIR" >&2; exit 3; }

# Soft-masked assembly so the CDS-masking / TE code paths are exercised.
SCAFFOLDS="${SCAFFOLDS:-$(find_one "${GENOME}_AssemblyScaffolds_Repeatmasked.fasta.gz")}"
[[ -z "$SCAFFOLDS" ]] && SCAFFOLDS="$(find_one "${GENOME}_AssemblyScaffolds.fasta.gz")"
[[ -n "$SCAFFOLDS" && -f "$SCAFFOLDS" ]] || { echo "ERROR: scaffolds not found under $RAW_DIR" >&2; exit 3; }

# Domains (IPR/Pfam) — prefer the GeneCatalog table so its protein IDs match the
# GeneCatalog GFF3 above.
DOMAINS="${DOMAINS:-$(find_one "${GENOME}_GeneCatalog_*_IPR.tab.gz")}"
[[ -z "$DOMAINS" ]] && DOMAINS="$(find_one "${GENOME}_*_IPR.tab.gz")"

# Quick sanity probe: the default -md mapper expects the 'HMMPfam' method label.
# Older JGI IPR tables sometimes use 'Pfam' instead — warn so a 0 TE factor
# isn't mistaken for parity.
if [[ -n "$DOMAINS" ]] && command -v zcat >/dev/null 2>&1; then
  if ! zcat "$DOMAINS" 2>/dev/null | head -n 200 | grep -q "HMMPfam"; then
    echo "NOTE: '$DOMAINS' has no 'HMMPfam' label in its first lines." >&2
    echo "      If the TE factor comes out 0, pass a custom -md regex (peek with: zcat <file> | head)." >&2
  fi
fi

# Shared argument list — byte-for-byte identical for both tools.
COMMON_ARGS=( -ig "$GFF3" -is "$SCAFFOLDS" -ib "$BUSCO_STRING" )
if [[ -n "$DOMAINS" ]]; then
  COMMON_ARGS+=( -id "$DOMAINS" )
else
  echo "WARNING: no ${GENOME}_*_IPR.tab.gz found — TE factor and domain counts will be ~0 in BOTH tools." >&2
fi

echo "==> Inputs:"
printf '    %s\n' "GFF3:      $GFF3" "SCAFFOLDS: $SCAFFOLDS" "DOMAINS:   ${DOMAINS:-<none>}" "BUSCO:     $BUSCO_STRING"

JAVA_OUT="$OUT_DIR/${GENOME}_java.tsv"
PY_OUT="$OUT_DIR/${GENOME}_python.tsv"
REPORT="$OUT_DIR/${GENOME}_parity_report.tsv"

# ---------- Run JAVA ----------
echo "==> Running Java…"
j_start=$(date +%s.%N)
java -jar "$JAR" "${COMMON_ARGS[@]}" -o "$JAVA_OUT"
j_end=$(date +%s.%N)
J_SECS=$(awk "BEGIN{printf \"%.2f\", $j_end-$j_start}")

# ---------- Run PYTHON (same inputs; auto-run disabled so it matches Java) ----------
echo "==> Running Python…"
p_start=$(date +%s.%N)
$PY_ABAQS "${COMMON_ARGS[@]}" --no-busco-auto-run -o "$PY_OUT"
p_end=$(date +%s.%N)
P_SECS=$(awk "BEGIN{printf \"%.2f\", $p_end-$p_start}")

echo "==> Runtime — Java: ${J_SECS}s   Python: ${P_SECS}s"

# ---------- Compare: writes the combined report and sets the exit code ----------
status=0
python3 "$COMPARE" \
  --genome "$GENOME" \
  --java "$JAVA_OUT" --python "$PY_OUT" \
  --java-secs "$J_SECS" --python-secs "$P_SECS" \
  --tol "$TOL" \
  --report "$REPORT" || status=$?

echo
echo "==> Combined report written to: $REPORT"
exit $status
