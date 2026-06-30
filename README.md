# ABAQS (Python)

**Accuracy-Based Annotation Quality Score** — a Python 3 port of the LBL/JGI MycoCosm
ABAQS tool. ABAQS scores the quality of a eukaryotic genome annotation from its GFF3,
assembly, protein domains, and BUSCO completeness, producing seven factors and a single
composite score.

This port is validated to **exact numerical parity** with the original Java implementation
(see [Validation](#validation)) and adds first-class **compleasm** support with automatic
**BUSCO** fallback.

---

## Installation

ABAQS itself is pure-Python (stdlib + biopython). The optional BUSCO-completeness step uses
**compleasm** (which needs **miniprot**); BUSCO can be used as a fallback. The reproducible
path is a conda/mamba environment:

```bash
# 1. Environment with the native tools ABAQS shells out to
mamba create -n abaqs -c conda-forge -c bioconda python=3.12 miniprot biopython pandas pip
conda activate abaqs

# 2. ABAQS itself
pip install git+https://github.com/iPsychonaut/abaqs-py.git

# 3. compleasm — installed as its self-contained script.
#    (It is not on PyPI, and the bioconda recipe's sepp/pasta/dendropy pins do not solve;
#     compleasm's only real runtime dependency is miniprot, already installed above.)
curl -fsSL https://raw.githubusercontent.com/huangnengCSU/compleasm/v0.2.6/compleasm.py \
     -o "$CONDA_PREFIX/bin/compleasm" && chmod +x "$CONDA_PREFIX/bin/compleasm"

# 4. (optional) BUSCO as a fallback completeness engine
# mamba install -n abaqs -c conda-forge -c bioconda busco
```

> The `pip install` requires the corrected `build-backend` in this repo
> (`setuptools.build_meta`). Installing an older revision will fail with
> `Cannot import 'setuptools.backends.legacy'`.

A one-shot setup script that also builds the Java reference and the differential-test harness
lives in [`tests/parity/setup_env.sh`](tests/parity/setup_env.sh).

---

## Usage

```bash
abaqs -ig annotation.gff3.gz \
      -is assembly.fasta.gz \
      -id domains_IPR.tab.gz \
      -o results.tsv
```

Only `-ig` (GFF3) is strictly required. If proteins are not supplied (`-ip`) they are
translated from the CDS features; if a domains table is not supplied (`-id`) the TE factor
degenerates. Run `abaqs --help` for the full flag list (drop-in compatible with the Java CLI).

### Output

A 13-line tab-separated report: six counts plus the seven factors and the composite score.

```
Total records / genes / scaffolds / proteins / proteins with domains / unique domains
Protein lengths distrbution factor
Incomplete genes factor
Transposable elements factor
Isoforms factor
BUSCO duplicated factor
BUSCO complete factor
ABAQS score
```

---

## BUSCO completeness: compleasm by default, BUSCO as backup

ABAQS resolves the BUSCO `complete`/`duplicated` factors in this order; the first that
succeeds wins:

1. **`-ib "C:..%[S:..%,D:..%],...,n:758"`** — an explicit BUSCO summary string.
2. **`-ibf summary.txt`** — an explicit compleasm `summary.txt` or BUSCO `short_summary*.txt`
   (format auto-detected).
3. **Auto-discovery** next to `--input-scaffolds-fasta` (EGAP layout):
   `<basename>_<lineage>_compleasm/summary.txt`, then
   `<basename>_<lineage>_busco/short_summary*.txt`.
4. **Auto-run** (default): run **compleasm** on the assembly; **if that fails, run BUSCO**.
5. If nothing yields data, ABAQS warns and assumes `complete=1.0, duplicated=0.0`.

So out of the box ABAQS **runs compleasm by default and falls back to BUSCO**. Relevant flags:

| Flag | Default | Meaning |
|---|---|---|
| `--busco-lineage` | `eukaryota_odb10` | Lineage to discover/run (e.g. `fungi_odb10`) |
| `--busco-threads` | `4` | Threads passed to compleasm/BUSCO |
| `--no-busco-auto-run` | (off) | Disable discovery + auto-run; use only `-ib`/`-ibf` |
| `-ib` / `-ibf` | — | Provide BUSCO data directly (skips running anything) |

---

## Validation

The port is checked by a three-way differential harness in
[`tests/parity/`](tests/parity): **Java (fixed BUSCO)** vs **Python (fixed BUSCO)** vs
**Python (compleasm)** on identical inputs.

- **PARITY** — Java vs Python on a fixed BUSCO string must match exactly (display + numeric).
- **REGRESS** — Python fixed-BUSCO vs Python compleasm: the six counts and the four
  non-BUSCO factors must be identical; only the BUSCO factors and the score may move.

Result on *Alternaria alternata* `Altalt1` (13,086 genes; see [Test data](#test-data)):

```
metric                                  java     py-busco   py-compleasm    B vs C   PARITY  REGRESS
Total records                          13111        13111          13111        +0   OK      OK
Total genes                            13086        13086          13086        +0   OK      OK
Total scaffolds                           26           26             26        +0   OK      OK
Total proteins                         13086        13086          13086        +0   OK      OK
Total proteins with domains             8702         8702           8702        +0   OK      OK
Total unique domains                    3579         3579           3579        +0   OK      OK
Protein lengths distrbution factor    0.9543       0.9543         0.9543   +0.0000   OK      OK
Incomplete genes factor               0.9745       0.9745         0.9745   +0.0000   OK      OK
Transposable elements factor          0.0511       0.0511         0.0511   +0.0000   OK      OK
Isoforms factor                       0.0000       0.0000         0.0000   +0.0000   OK      OK
BUSCO duplicated factor               0.0040       0.0040         0.0018   -0.0022   OK      INFO
BUSCO complete factor                 0.9930       0.9930         0.9911   -0.0019   OK      INFO
ABAQS score                           0.9353       0.9353         0.9354   +0.0001   OK      INFO

PARITY  (java vs py-busco):      exact (tier1)
REGRESS (py-busco vs compleasm): non-BUSCO identical
VERDICT: PASS — port reproduces Java; compleasm perturbs only BUSCO outputs
```

(The fixed BUSCO string is the `fungi_odb10` placeholder `n:758`; compleasm resolved
`fungi_odb12`, N=1122 — hence the expected, informational shift in the BUSCO rows only.)

To reproduce:

```bash
bash tests/parity/setup_env.sh                              # builds the Java jar + env
RAW_DIR=/path/to/jgi/files bash tests/parity/run_all_tests.sh
```

### Test data

Validation uses the public JGI MycoCosm genome **_Alternaria alternata_ MPI-PUGE-AT-0064 v1.0**
(portal **`Altalt1`**, project AP-1103619). The files are **not redistributed here** (JGI data
requires a (free) account and is subject to the JGI Data Use Policy) — download them from the
MycoCosm `Altalt1` portal:

| File | MycoCosm location | Role |
|---|---|---|
| `Altalt1_GeneCatalog_20170627.gff3.gz` | Annotation ▸ Filtered Models (best) ▸ Genes | GFF3 (`-ig`) |
| `Altalt1_AssemblyScaffolds_Repeatmasked.fasta.gz` | Assembly ▸ Genome Assembly (masked) | soft-masked scaffolds (`-is`) |
| `Altalt1_GeneCatalog_proteins_20170627_IPR.tab.gz` | Annotation ▸ Filtered Models (best) ▸ Functional Annotations ▸ InterPro | Pfam/InterPro domains (`-id`) |
| `Altalt1_GeneCatalog_proteins_20170627.aa.fasta.gz` *(optional)* | Annotation ▸ Filtered Models (best) ▸ Proteins | pre-translated proteins (`-ip`) |

Keep all inputs from the same **Filtered Models (best)** set so protein IDs line up across the
GFF3 and the domain table.

---

## License

See [LICENSE](LICENSE). ABAQS is a port of the LBL/JGI MycoCosm tool; the scoring methodology
is the original authors'.
