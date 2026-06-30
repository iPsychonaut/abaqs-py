#!/usr/bin/env bash
#
# ONE test that runs all three and writes a single combined report:
#   1. Java        on (GFF3, scaffolds, domains) + fixed BUSCO string
#   2. Python port on the SAME inputs            + fixed BUSCO string
#   3. Python port on the SAME inputs            + REAL compleasm BUSCO
#
# Checks (see compare_three.py):
#   PARITY   = Java vs Python(fixed)      -> the port must reproduce Java.
#   REGRESS  = Python(fixed) vs compleasm -> only BUSCO factors + score may move.
#
# Prereqs: run setup_env.sh once, then `conda activate abaqs` (the env holds
# both the port and compleasm). Put the JGI files in $RAW_DIR.
#
set -euo pipefail

# ---------- Configuration (override via env) ----------
ABAQS_JAVA_DIR="${ABAQS_JAVA_DIR:-/mnt/c/Users/theda/Desktop/abaqs-java}"
JAR="${JAR:-$ABAQS_JAVA_DIR/target/abaqs-jar-with-dependencies.jar}"
PY_ABAQS="${PY_ABAQS:-abaqs}"
ABAQS_ENV="${ABAQS_ENV:-abaqs}"
# Work even without `conda activate abaqs`: fall back to the env's abaqs binary
# (its pip shebang pins the env python, so the port imports correctly).
if ! command -v "$PY_ABAQS" >/dev/null 2>&1; then
  cand="$(conda info --base 2>/dev/null)/envs/$ABAQS_ENV/bin/abaqs"
  [ -x "$cand" ] && PY_ABAQS="$cand"
fi
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
COMPARE="${COMPARE:-$SCRIPT_DIR/compare_three.py}"

RAW_DIR="${RAW_DIR:-$HOME/abaqs-parity/raw}"
OUT_DIR="${OUT_DIR:-$HOME/abaqs-parity/out}"
GENOME="${GENOME:-Altalt1}"
TOL="${TOL:-0.0001}"
LINEAGE="${LINEAGE:-fungi_odb10}"               # matches the placeholder n:758
THREADS="${BUSCO_THREADS:-4}"
BUSCO_STRING="${BUSCO_STRING:-C:99.3%[S:98.9%,D:0.4%],F:0.3%,M:0.4%,n:758}"
COMPLEASM_SUMMARY="${COMPLEASM_SUMMARY:-}"      # precomputed summary.txt (optional)
COMPLEASM_LIB="${COMPLEASM_LIB:-}"              # pre-downloaded lineage library (-L, optional)

mkdir -p "$OUT_DIR"; WORK="$OUT_DIR/work"; mkdir -p "$WORK"
[[ -f "$JAR" ]] || { echo "ERROR: jar not found: $JAR  (run setup_env.sh)" >&2; exit 2; }

now()  { date +%s.%N; }
secs() { awk "BEGIN{printf \"%.2f\", $2-$1}"; }

# ---------- Resolve inputs ----------
find_one() { find "$RAW_DIR" -type f -iname "$1" 2>/dev/null | sort | head -n1; }

GFF3="${GFF3:-$(find_one "${GENOME}_GeneCatalog_*.gff3.gz")}"
[[ -z "$GFF3" ]] && GFF3="$(find_one "${GENOME}_GeneCatalog_*.gff3")"
[[ -n "$GFF3" && -f "$GFF3" ]] || { echo "ERROR: GFF3 not found under $RAW_DIR" >&2; exit 3; }

SCAFFOLDS="${SCAFFOLDS:-$(find_one "${GENOME}_AssemblyScaffolds_Repeatmasked.fasta.gz")}"
[[ -z "$SCAFFOLDS" ]] && SCAFFOLDS="$(find_one "${GENOME}_AssemblyScaffolds.fasta.gz")"
[[ -n "$SCAFFOLDS" && -f "$SCAFFOLDS" ]] || { echo "ERROR: scaffolds not found under $RAW_DIR" >&2; exit 3; }

DOMAINS="${DOMAINS:-$(find_one "${GENOME}_GeneCatalog_*_IPR.tab.gz")}"
[[ -z "$DOMAINS" ]] && DOMAINS="$(find_one "${GENOME}_*_IPR.tab.gz")"

# Shared scoring args (identical across all three runs, sans BUSCO source).
SCORE_ARGS=( -ig "$GFF3" -is "$SCAFFOLDS" )
[[ -n "$DOMAINS" ]] && SCORE_ARGS+=( -id "$DOMAINS" ) \
  || echo "WARNING: no IPR domains file — TE factor/domain counts ~0 in all runs." >&2

echo "==> Inputs:"
printf '    %s\n' "GFF3:      $GFF3" "SCAFFOLDS: $SCAFFOLDS" "DOMAINS:   ${DOMAINS:-<none>}"

JAVA_OUT="$OUT_DIR/${GENOME}_java.tsv"
PB_OUT="$OUT_DIR/${GENOME}_python_busco.tsv"
PC_OUT="$OUT_DIR/${GENOME}_python_compleasm.tsv"
REPORT="$OUT_DIR/${GENOME}_three_way_report.tsv"

# ---------- 1) Java (fixed BUSCO) ----------
echo "==> [1/3] Java (fixed BUSCO)…"
t0=$(now); java -jar "$JAR" "${SCORE_ARGS[@]}" -ib "$BUSCO_STRING" -o "$JAVA_OUT"; t1=$(now)
J_SECS=$(secs "$t0" "$t1")

# ---------- 2) Python (fixed BUSCO) ----------
echo "==> [2/3] Python (fixed BUSCO)…"
t0=$(now)
$PY_ABAQS "${SCORE_ARGS[@]}" -ib "$BUSCO_STRING" --no-busco-auto-run -o "$PB_OUT"
t1=$(now); PB_SECS=$(secs "$t0" "$t1")

# ---------- obtain a compleasm summary.txt (reuse if already computed) ----------
CDIR="$OUT_DIR/${GENOME}_${LINEAGE}_compleasm"
if [[ -n "$COMPLEASM_SUMMARY" && -f "$COMPLEASM_SUMMARY" ]]; then
  SUMMARY="$COMPLEASM_SUMMARY"
  echo "==> Using precomputed compleasm summary: $SUMMARY"
elif [[ -s "$CDIR/summary.txt" ]]; then
  SUMMARY="$CDIR/summary.txt"
  echo "==> Reusing existing compleasm summary: $SUMMARY"
else
  # Locate the compleasm binary: explicit > PATH > any existing conda env.
  if [[ -z "${COMPLEASM_BIN:-}" ]]; then
    COMPLEASM_BIN="$(command -v compleasm 2>/dev/null || true)"
    if [[ -z "$COMPLEASM_BIN" ]] && command -v conda >/dev/null 2>&1; then
      cbase="$(conda info --base 2>/dev/null)"
      COMPLEASM_BIN="$(ls "$cbase"/bin/compleasm "$cbase"/envs/*/bin/compleasm 2>/dev/null | head -n1)"
    fi
  fi
  if [[ -z "$COMPLEASM_BIN" ]]; then
    echo "ERROR: compleasm not found. Set COMPLEASM_BIN=/path, activate an env" >&2
    echo "       that has it, or pass COMPLEASM_SUMMARY=/path/to/summary.txt." >&2
    exit 4
  fi
  GENOME_FA="$WORK/${GENOME}_genome.fasta"        # compleasm/miniprot want plain FASTA
  if [[ "$SCAFFOLDS" == *.gz ]]; then
    zcat "$SCAFFOLDS" > "$GENOME_FA"
  else
    cp -f "$SCAFFOLDS" "$GENOME_FA"
  fi
  # Run compleasm with the env's OWN python (it imports pandas) and the env bin
  # on PATH (so it finds miniprot); explicit python bypasses the script shebang,
  # which would otherwise pick the system python.
  CL_DIR="$(dirname "$COMPLEASM_BIN")"
  CL_PY="$CL_DIR/python"
  [[ -x "$CL_PY" ]] || CL_PY="$(command -v python3)"
  cl=(run -a "$GENOME_FA" -o "$CDIR" -l "$LINEAGE" -t "$THREADS")
  [[ -n "$COMPLEASM_LIB" ]] && cl+=(-L "$COMPLEASM_LIB")
  echo "==> Running compleasm ($COMPLEASM_BIN, lineage=$LINEAGE, threads=$THREADS)…"
  PATH="$CL_DIR:$PATH" "$CL_PY" "$COMPLEASM_BIN" "${cl[@]}"
  SUMMARY="$CDIR/summary.txt"
fi
[[ -f "$SUMMARY" ]] || { echo "ERROR: compleasm summary not found: $SUMMARY" >&2; exit 4; }

# ---------- 3) Python (compleasm BUSCO) ----------
echo "==> [3/3] Python (compleasm BUSCO)…"
t0=$(now); $PY_ABAQS "${SCORE_ARGS[@]}" -ibf "$SUMMARY" -o "$PC_OUT"; t1=$(now)
PC_SECS=$(secs "$t0" "$t1")

echo "==> Runtime — java:${J_SECS}s  py-busco:${PB_SECS}s  py-compleasm:${PC_SECS}s"

# ---------- compare all three ----------
status=0
python3 "$COMPARE" \
  --genome "$GENOME" \
  --java "$JAVA_OUT" --py-busco "$PB_OUT" --py-compleasm "$PC_OUT" \
  --secs-java "$J_SECS" --secs-py-busco "$PB_SECS" --secs-py-compleasm "$PC_SECS" \
  --tol "$TOL" --report "$REPORT" || status=$?

echo
echo "==> Three-way report: $REPORT"
exit $status
