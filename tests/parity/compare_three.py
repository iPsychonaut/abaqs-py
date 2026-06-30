#!/usr/bin/env python3
"""Three-way ABAQS comparison in one report:

  A = Java        (fixed BUSCO string)
  B = Python      (fixed BUSCO string)   -- same inputs as A
  C = Python      (compleasm BUSCO)      -- same inputs, BUSCO from compleasm

Two independent checks:
  * PARITY   (A vs B): the port must reproduce Java on every field.
       - tier1 = exact formatted match; tier2 = ints equal, floats within --tol.
  * REGRESS  (B vs C): switching the BUSCO source must change ONLY the 2 BUSCO
       factors + ABAQS score. The 6 counts + 4 non-BUSCO factors must match;
       the BUSCO-driven rows are informational (reported with deltas, never fail).

Exit 0 iff PARITY holds (tier2; tier1 with --strict) AND REGRESS holds over the
non-BUSCO metrics. Else exit 1.
"""
from __future__ import annotations

import argparse
import sys

INT_KEYS = [
    "Total records", "Total genes", "Total scaffolds",
    "Total proteins", "Total proteins with domains", "Total unique domains",
]
FLOAT_KEYS = [
    "Protein lengths distrbution factor",   # NB: 'distrbution' typo preserved
    "Incomplete genes factor",
    "Transposable elements factor",
    "Isoforms factor",
    "BUSCO duplicated factor",
    "BUSCO complete factor",
    "ABAQS score",
]
ORDER = INT_KEYS + FLOAT_KEYS
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


def _num_eq(key, x, y, tol):
    try:
        if key in INT_KEYS:
            return int(x) == int(y)
        return abs(float(x) - float(y)) <= tol
    except (ValueError, TypeError):
        return False


def _delta(key, x, y):
    try:
        if key in INT_KEYS:
            return f"{int(y) - int(x):+d}"
        return f"{float(y) - float(x):+.4f}"
    except (ValueError, TypeError):
        return ""


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--genome", default="")
    ap.add_argument("--java", required=True)
    ap.add_argument("--py-busco", required=True)
    ap.add_argument("--py-compleasm", required=True)
    ap.add_argument("--report", required=True)
    ap.add_argument("--tol", type=float, default=1e-4)
    ap.add_argument("--secs-java", default="?")
    ap.add_argument("--secs-py-busco", default="?")
    ap.add_argument("--secs-py-compleasm", default="?")
    ap.add_argument("--strict", action="store_true",
                    help="PARITY requires exact (tier1) match for exit 0")
    a = ap.parse_args()

    A, B, C = parse(a.java), parse(a.py_busco), parse(a.py_compleasm)

    rows = []
    parity_t1 = parity_t2 = regress_ok = True
    for key in ORDER:
        av, bv, cv = A.get(key), B.get(key), C.get(key)
        info = key in BUSCO_INFO_KEYS

        # PARITY: Java (A) vs Python fixed-BUSCO (B)
        if av is None or bv is None:
            p_stat, p_t1, p_t2 = "MISSING", False, False
        else:
            p_t1 = (av == bv)
            p_t2 = _num_eq(key, av, bv, a.tol)
            p_stat = "OK" if p_t1 else ("~tol" if p_t2 else "DIFF")
        parity_t1 = parity_t1 and p_t1
        parity_t2 = parity_t2 and p_t2

        # REGRESS: Python fixed (B) vs Python compleasm (C)
        if bv is None or cv is None:
            r_stat, r_ok = "MISSING", False
        elif info:
            r_stat, r_ok = "INFO", True            # expected to move
        else:
            r_ok = _num_eq(key, bv, cv, a.tol)
            r_stat = "OK" if r_ok else "DIFF"
        if not info:
            regress_ok = regress_ok and r_ok

        rows.append((key, av, bv, cv, _delta(key, bv, cv), p_stat, r_stat))

    parity_pass = parity_t1 if a.strict else parity_t2
    overall = parity_pass and regress_ok
    code = 0 if overall else 1

    parity_word = ("exact (tier1)" if parity_t1
                   else (f"within tol={a.tol} (tier2)" if parity_t2 else "FAILED"))
    regress_word = "non-BUSCO identical" if regress_ok else "FAILED"

    # ---- combined report ----
    w = max(len(k) for k in ORDER)
    out = []
    out.append(f"# ABAQS three-way report   genome={a.genome}   tol={a.tol}")
    out.append(f"# runtime  java={a.secs_java}s  py-busco={a.secs_py_busco}s  "
               f"py-compleasm={a.secs_py_compleasm}s")
    out.append("")
    for title, path in (("JAVA (fixed BUSCO)", a.java),
                        ("PYTHON (fixed BUSCO)", a.py_busco),
                        ("PYTHON (compleasm)", a.py_compleasm)):
        out.append(f"=== {title} ===")
        out.append(open(path, encoding="utf-8").read().rstrip("\n"))
        out.append("")
    out.append("=== COMPARISON ===")
    out.append(f"{'metric'.ljust(w)}\t{'java':>13}\t{'py-busco':>13}\t"
               f"{'py-compleasm':>13}\t{'Bvs C':>8}\tPARITY\tREGRESS")
    for key, av, bv, cv, dl, p_stat, r_stat in rows:
        out.append(f"{key.ljust(w)}\t{str(av):>13}\t{str(bv):>13}\t"
                   f"{str(cv):>13}\t{dl:>8}\t{p_stat}\t{r_stat}")
    out.append("")
    out.append(f"PARITY  (java vs py-busco):     {parity_word}")
    out.append(f"REGRESS (py-busco vs compleasm): {regress_word}")
    out.append(f"VERDICT: {'PASS' if overall else 'FAIL'} — "
               + ("port reproduces Java; compleasm perturbs only BUSCO outputs"
                  if overall else "see DIFF rows above"))

    text = "\n".join(out) + "\n"
    with open(a.report, "w", encoding="utf-8") as fh:
        fh.write(text)
    sys.stdout.write(text)
    return code


if __name__ == "__main__":
    sys.exit(main())
