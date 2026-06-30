# ABAQS (Python)

**Accuracy-Based Annotation Quality Score** — a Python 3 port/refactor of the original
**ABAQS** tool developed for the LBL/JGI MycoCosm annotation-quality workflow.

ABAQS scores the quality of a eukaryotic genome annotation from its GFF3, assembly,
protein-domain annotations, and BUSCO completeness, producing six summary counts,
seven quality factors, and one composite score.

This repository preserves the original ABAQS scoring logic while providing a Python
implementation, reproducibility tests, and **compleasm** support with automatic
**BUSCO** fallback. The Python implementation is validated to **exact numerical parity**
with the original Java implementation when the same fixed BUSCO summary is supplied
(see [Validation](#validation)).

## Relationship to the original ABAQS

The original ABAQS implementation is maintained at
[mycocosm-lbl/abaqs](https://github.com/mycocosm-lbl/abaqs). The original ABAQS software
and scoring methodology were developed by the original LBL/JGI MycoCosm authors. This
repository is a Python port/refactor intended to reproduce that behavior and make the
tool easier to install, test, and integrate into Python-centered annotation workflows.

This repository is not an official LBL/JGI, Lawrence Berkeley National Laboratory,
University of California, U.S. Department of Energy, or MycoCosm release. Those names
are used only to identify the provenance of the original ABAQS tool and should not be
interpreted as endorsement of this Python port.

## Table of Contents

1. [Relationship to the original ABAQS](#relationship-to-the-original-abaqs)
2. [Installation](#installation)
3. [Usage](#usage)
4. [Input requirements and Java compatibility](#input-requirements-and-java-compatibility)
5. [Interpreting the results](#interpreting-the-results)
6. [BUSCO completeness: compleasm by default, BUSCO as backup](#busco-completeness-compleasm-by-default-busco-as-backup)
7. [Validation](#validation)
8. [Test data](#test-data)
9. [Attribution and license](#attribution-and-license)

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
Protein lengths distribution factor
Incomplete genes factor
Transposable elements factor
Isoforms factor
BUSCO duplicated factor
BUSCO complete factor
ABAQS score

```

The output label `Protein lengths distribution factor` is intentionally preserved for
drop-in display compatibility with the original Java output.

---

## Input requirements

Basic ABAQS scoring is designed around four major inputs:

1. **Assembly FASTA (`-is`)**
   - The assembly FASTA should generally be soft-masked with a repeat-masking program.
   - Low-complexity repeats are ideally ignored during repeat masking, for example, by
     running RepeatMasker with `-nolow` when appropriate.

2. **Gene models in GFF3 format (`-ig`)**
   - This is the only strictly required input in `abaqs-py`.
   - For full Java-compatible scoring, the GFF3 should contain protein-coding gene models
     with stable protein identifiers.
   - The original Java implementation expects protein identifiers in the GFF3 attributes
     field as `proteinId` by default. Use `-mg` to specify an alternate mapper.

3. **BUSCO or compleasm completeness data (`-ib` or `-ibf`)**
   - BUSCO data can be supplied directly as a summary string with `-ib`.
   - A BUSCO `short_summary*.txt` or compleasm `summary.txt` file can be supplied with
     `-ibf`.
   - `abaqs-py` can also discover or run compleasm/BUSCO automatically when assembly
     inputs and lineage settings are available.

4. **Protein-domain data (`-id`)**
   - Domain data are used for the transposable-element and domain-loading components.
   - The original ABAQS default expects a tab-separated domain file where the protein ID
     can be mapped from the domain records and the domain source contains `HMMPfam`.
   - If using InterProScan TSV output where the source field contains `Pfam`, use an
     appropriate `-md` mapper.

Optional protein FASTA input can be supplied with `-ip`. If proteins are not supplied,
`abaqs-py` translates proteins from CDS features when possible.

abaqs-py aims to remain drop-in compatible with the original Java CLI where practical,
while also adding Python-specific testing and compleasm support.

---

## Interpreting the results

After a successful run, ABAQS reports six counts, seven factors, and the final score.

| Output row | Meaning |
|---|---|
| `Total records` | Total number of features in the GFF3 file. |
| `Total genes` | Count of protein-coding genes detected in the GFF3. |
| `Total scaffolds` | Number of sequence records in the assembly FASTA. |
| `Total proteins` | Number of proteins used for scoring. Ideally, this matches the protein-coding gene count. |
| `Total proteins with domains` | Number of proteins with at least one mapped Pfam/InterPro domain. Very low values may indicate incomplete domain annotation or a domain-file parsing issue. |
| `Total unique domains` | Number of unique domain identifiers in the supplied domain table. |
| `Protein lengths distribution factor` | Protein-length distribution factor. Values closer to 1 are better; low values reduce the final ABAQS score. |
| `Incomplete genes factor` | Proportion of genes with expected start/stop completeness. Low values reduce the final score. |
| `Transposable elements factor` | Proportion of retained predicted proteins classified as transposable-element-associated according to domain and masking criteria. |
| `Isoforms factor` | Fraction of isoform-like overlapping coding sequences. The original ABAQS scoring expects this to be low. |
| `BUSCO duplicated factor` | BUSCO-derived duplication estimate. |
| `BUSCO complete factor` | BUSCO- or compleasm-derived completeness estimate. |
| `ABAQS score` | Final composite ABAQS score on a 0–1 scale. |

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

## Attribution and license

abaqs-py is a Python port/refactor of the original ABAQS tool from the LBL/JGI MycoCosm
ecosystem. The original ABAQS scoring methodology and Java reference implementation are
the work of the original ABAQS authors.

The original ABAQS repository is:

- https://github.com/mycocosm-lbl/abaqs

The original ABAQS license permits redistribution and modification, but requires that
redistributions retain the original copyright notice, license conditions, disclaimer, and
no-endorsement clause. The original license also contains an additional enhancement-grant
clause. See [LICENSE](LICENSE) for the retained original ABAQS license notice and the
additional abaqs-py copyright notice.

abaqs-py adds:

- a Python 3 implementation,
- Python packaging,
- parity tests against the original Java implementation,
- compleasm support,
- BUSCO fallback/discovery logic,
- repository-specific test harnesses and documentation.

This repository is not an official LBL/JGI, Lawrence Berkeley National Laboratory,
University of California, U.S. Department of Energy, or MycoCosm release. Those names
may not be used to endorse or promote this project without specific prior written
permission.
