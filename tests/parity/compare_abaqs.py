#!/usr/bin/env python3
"""Two-tier differential comparison for ABAQS output files.

Default mode (Java vs Python parity):
  * Tier 1 (strict): the formatted output values must be byte-identical.
  * Tier 2 (tolerant): integer counts equal; float factors within --tol.

--expect-busco-diff mode (fixed-BUSCO vs compleasm, both Python):
  * the 6 counts + 4 non-BUSCO factors MUST still match (regression guard);
  * the 2 BUSCO factors + ABAQS score are EXPECTED to differ and are reported
    as informational rows (with deltas) that never cause a failure.

Writes a single combined report with BOTH inputs' raw output plus the table.

Exit code:
  0  -> the (non-informational) metrics pass Tier 2.  With --strict, requires
        Tier 1 (exact) over those metrics.
  1  -> a compared metric diverges (or is missing).
"""
from __future__ import annotations

import argparse
import sys

INT_KEYS = [
    "Total records",
    "Total genes",
    "Total scaffolds",
    "Total proteins",
    "Total proteins with domains",
    "Total unique domains",
]
FLOAT_KEYS = [
    "Protein lengths distrbution factor",   # NB: 'distrbution' typo preserved in both tools
    "Incomplete genes factor",
    "Transposable elements factor",
    "Isoforms factor",
    "BUSCO duplicated factor",
    "BUSCO complete factor",
    "ABAQS score",
]
ORDER = INT_KEYS + FLOAT_KEYS

# Outputs that depend on BUSCO — expected to move when the BUSCO source changes.
BUSCO_INFO_KEYS = {"BUSCO duplicated factor", "BUSCO complete factor", "ABAQS score"}


def parse(path: str) -> dict:
    d = {}
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            if "\t" not in line:
                continue
            key, val = line.rstrip("\n").split("\t", 1)
            d[key.rstrip(":").strip()] = val.strip()
    return d


def _delta(key: str, av: str, bv: str) -> str:
    try:
        if key in INT_KEYS:
            return f"{int(bv) - int(av):+d}"
        return f"{float(bv) - float(av):+.4f}"
    except ValueError:
        return ""


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--genome", default="")
    ap.add_argument("--java", required=True, help="file A")
    ap.add_argument("--python", required=True, help="file B")
    ap.add_argument("--label-a", default="java")
    ap.add_argument("--label-b", default="python")
    ap.add_argument("--report", required=True)
    ap.add_argument("--tol", type=float, default=1e-4)
    ap.add_argument("--java-secs", default="?")
    ap.add_argument("--python-secs", default="?")
    ap.add_argument("--strict", action="store_true",
                    help="require exact (Tier 1) parity for exit 0")
    ap.add_argument("--expect-busco-diff", action="store_true",
                    help="treat BUSCO factors + ABAQS score as informational")
    a = ap.parse_args()

    A = parse(a.java)
    B = parse(a.python)

    rows = []
    tier1_ok = True
    tier2_ok = True
    for key in ORDER:
        av, bv = A.get(key), B.get(key)
        info = a.expect_busco_diff and key in BUSCO_INFO_KEYS
        if av is None or bv is None:
            rows.append((key, av or "<missing>", bv or "<missing>", "MISSING", "MISSING", "", info))
            if not info:
                tier1_ok = tier2_ok = False
            continue
        t1 = (av == bv)
        if key in INT_KEYS:
            try:
                t2 = int(av) == int(bv)
            except ValueError:
                t2 = False
        else:
            try:
                t2 = abs(float(av) - float(bv)) <= a.tol
            except ValueError:
                t2 = False
        if not info:
            tier1_ok = tier1_ok and t1
            tier2_ok = tier2_ok and t2
        s1 = "INFO" if info else ("OK" if t1 else "DIFF")
        s2 = "INFO" if info else ("OK" if t2 else "DIFF")
        rows.append((key, av, bv, s1, s2, _delta(key, av, bv), info))

    # ---- verdict ----
    if a.expect_busco_diff:
        if tier2_ok:
            verdict = ("PASS — non-BUSCO metrics identical; BUSCO factors + ABAQS "
                       "score differ as expected (compleasm vs fixed string)")
            code = 0
        else:
            verdict = f"FAIL — a non-BUSCO metric diverged beyond tol={a.tol}"
            code = 1
    elif tier1_ok:
        verdict, code = "PASS — exact parity (Tier 1)", 0
    elif tier2_ok:
        verdict = (f"PASS — numeric parity within tol={a.tol} (Tier 2); "
                   "Tier 1 differences are last-digit rounding only")
        code = 1 if a.strict else 0
    else:
        verdict, code = f"FAIL — numeric divergence beyond tol={a.tol}", 1

    # ---- combined report ----
    w = max(len(k) for k in ORDER)
    la, lb = a.label_a, a.label_b
    out = []
    out.append(f"# ABAQS comparison report   genome={a.genome}")
    out.append(f"# {la}={a.java_secs}s  {lb}={a.python_secs}s   tol={a.tol}"
               + ("   mode=expect-busco-diff" if a.expect_busco_diff else ""))
    out.append("")
    out.append(f"=== {la.upper()} OUTPUT ===")
    out.append(open(a.java, encoding="utf-8").read().rstrip("\n"))
    out.append("")
    out.append(f"=== {lb.upper()} OUTPUT ===")
    out.append(open(a.python, encoding="utf-8").read().rstrip("\n"))
    out.append("")
    out.append("=== COMPARISON ===")
    out.append(f"{'metric'.ljust(w)}\t{la:>14}\t{lb:>14}\t{'delta':>9}\tT1\tT2")
    for key, av, bv, s1, s2, dl, _info in rows:
        out.append(f"{key.ljust(w)}\t{av:>14}\t{bv:>14}\t{dl:>9}\t{s1}\t{s2}")
    out.append("")
    out.append(f"VERDICT: {verdict}")

    text = "\n".join(out) + "\n"
    with open(a.report, "w", encoding="utf-8") as fh:
        fh.write(text)
    sys.stdout.write(text)
    return code


if __name__ == "__main__":
    sys.exit(main())
