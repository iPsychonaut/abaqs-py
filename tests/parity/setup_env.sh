#!/usr/bin/env bash
#
# One-time setup for ABAQS parity testing — ONE conda env holds everything.
#   * apt:   JDK + Maven, builds the Java jar.
#   * mamba: a single env ('abaqs' by default) with Python 3.12, compleasm,
#            biopython, and the editable-installed Python port.
#
# Usage:
#   bash setup_env.sh
#   ABAQS_ENV=myenv bash setup_env.sh     # install into a different/existing env
#
set -euo pipefail

ABAQS_JAVA_DIR="${ABAQS_JAVA_DIR:-/mnt/c/Users/theda/Desktop/abaqs-java}"
ABAQS_PY_DIR="${ABAQS_PY_DIR:-/mnt/c/Users/theda/Desktop/abaqs-py}"
ABAQS_ENV="${ABAQS_ENV:-abaqs}"

echo "==> Installing the JDK + Maven (apt) to build the Java jar…"
sudo apt-get update
sudo apt-get install -y openjdk-21-jdk maven tar gzip findutils gawk curl

echo
echo "==> Building the Java jar…"
( cd "$ABAQS_JAVA_DIR" && mvn -q clean package )
JAR="$ABAQS_JAVA_DIR/target/abaqs-jar-with-dependencies.jar"
[ -f "$JAR" ] && echo "    OK: $JAR" || { echo "    ERROR: jar not built" >&2; exit 1; }

echo
command -v mamba >/dev/null 2>&1 || { echo "ERROR: mamba not on PATH." >&2; exit 1; }
# bioconda's `compleasm` recipe drags in sepp->pasta->dendropy, whose pins are
# self-contradictory and won't solve; and compleasm isn't on PyPI. But compleasm
# is a single self-contained script whose only external tool is miniprot. So we
# get miniprot + pandas from conda (trivial solve) and drop compleasm.py into the
# env's bin below. --override-channels avoids the Anaconda 'defaults' TOS channel.
CHANNELS=(--override-channels -c conda-forge -c bioconda)
if conda env list | awk '{print $1}' | grep -qx "$ABAQS_ENV"; then
  echo "==> Env '$ABAQS_ENV' exists — adding miniprot + pandas + biopython…"
  mamba install -y -n "$ABAQS_ENV" "${CHANNELS[@]}" miniprot pandas biopython
else
  echo "==> Creating env '$ABAQS_ENV' (python 3.12 + miniprot + pandas + biopython)…"
  mamba create -y -n "$ABAQS_ENV" "${CHANNELS[@]}" python=3.12 miniprot pandas biopython pip
fi

echo
echo "==> Installing the Python port (editable) into '$ABAQS_ENV'…"
conda run -n "$ABAQS_ENV" python -m pip install -e "$ABAQS_PY_DIR"

echo
echo "==> Installing the standalone compleasm script into the env's bin…"
CL_VER="${COMPLEASM_VERSION:-0.2.6}"
ENV_BIN="$(conda info --base)/envs/$ABAQS_ENV/bin"
curl -fsSL "https://raw.githubusercontent.com/huangnengCSU/compleasm/v${CL_VER}/compleasm.py" \
     -o "$ENV_BIN/compleasm" \
  || curl -fsSL "https://raw.githubusercontent.com/huangnengCSU/compleasm/main/compleasm.py" \
     -o "$ENV_BIN/compleasm"
chmod +x "$ENV_BIN/compleasm"

echo
echo "==> Verifying:"
conda run -n "$ABAQS_ENV" python --version
conda run -n "$ABAQS_ENV" abaqs --version
conda run -n "$ABAQS_ENV" compleasm --version 2>/dev/null || true

echo
echo "==> Setup complete. Everything (abaqs + compleasm) lives in env '$ABAQS_ENV'."
echo "    Run all three tests with:"
echo "        conda activate $ABAQS_ENV"
echo "        RAW_DIR=/mnt/c/Users/theda/Downloads bash run_all_tests.sh"
