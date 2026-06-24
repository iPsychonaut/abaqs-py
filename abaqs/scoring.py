"""Core ABAQS scoring pipeline and BUSCO summary parsing.

Mirrors:
  org.mycocosm.abaqs.main.ABAQS#processInput  (the orchestration + factors)
  org.mycocosm.abaqs.main.BuscoData           (BUSCO summary parsing)
  org.mycocosm.abaqs.main.GeneRecord          (per-mRNA record dataclass)
  org.mycocosm.abaqs.main.MapperDefinitionsFactory (mapper helpers)
  org.mycocosm.framework.text.PatternTransformer  (mapper pattern syntax)
  org.mycocosm.abaqs.main.PfamDomain          (collapsed to plain strings)

Public entry point:
    compute_abaqs(...)  →  ABAQSResult

Intentional divergences from Java (Python is more correct or more compatible):
  - A1: Pfam IDs are compared case-insensitively (uppercased before set lookup).
        Java's PfamDomain.equals is case-sensitive, so a Pfam list of "PF00078"
        will not match a domains file emitting "pf00078" in Java.
  - A2: mRNAs with no associated protein are silently skipped from the
        incomplete-genes factor (Java would NullPointerException).
  - A3/A4: GFF3 mapper accepts both bare attribute names and Java's
        ``column:attr:pattern->template`` form; FASTA mapper supports named
        groups and the ``<<CI,DA`` flag-suffix syntax. Java named-group syntax
        ``(?<name>…)`` is auto-translated to Python's ``(?P<name>…)``.
  - B1: BUSCO string parser respects raw counts when ``%`` is absent (Java
        unconditionally divides by 100 due to a regex-group nullability bug).

Python-only extensions (no Java equivalent):
  - compleasm support: ``BuscoData.of_compleasm_file`` parses compleasm's
    ``summary.txt``; ABAQS uses ``complete = S + D``.
  - Auto-discovery: when neither --busco-data nor --busco-data-file is
    provided, ABAQS searches for ``<basename>_<lineage>_compleasm/summary.txt``
    and ``<basename>_<lineage>_busco/short_summary*.txt`` beside the
    --input-scaffolds-fasta (EGAP layout). If nothing is found, ABAQS tries to
    run ``compleasm run``, then ``busco`` as a subprocess; on failure it warns
    and falls back to ``complete=1.0, duplicated=0.0`` (the 'no BUSCO data'
    behavior).

Algorithmic notes carried over from Java:
  - C3 (protein-length distribution): organism distribution is keyed by bin
        INDEX (0, 1, 2, …), while the reference distribution is keyed by the
        integer in column 1 of ``reference-proteins-length-distribution.tsv.gz``.
        The factor formula treats both as the same scalar (``length + 0.5``).
        Bundled reference data must therefore be in units consistent with
        ``--protein-length-binning`` (default 5 → keys are protein-length / 5).
"""

from __future__ import annotations

import gzip
import logging
import math
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Optional, Set

from abaqs.features import Feature, FeaturePair, FeatureTrack
from abaqs.sequence import (
    GFF3Data, Gff3Record, Gff3RecordCategory, Gff3Type,
    GeneCode,
    bin_collection_fixed_bin_size,
    load_default_gene_codes, load_fasta, load_gene_codes_file,
    open_gzipped_or_plain, parse_gff3_file,
)

log = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Defaults (mirror Java ABAQS constants)
# ---------------------------------------------------------------------------

DEFAULT_GFF3_PROTEIN_ID_ATTR         = 'proteinId'
DEFAULT_ISOFORMS_MIN_OVERLAP         = 0.25
DEFAULT_MASKER_FUNCTION              = 'TO_LOWER_CASE'
DEFAULT_NO_DOMAINS_CDS_MASKED_CUTOFF = 0.20
DEFAULT_SUSPECTED_DOMAINS_MASKED_CUTOFF = float('nan')  # NaN = always TE
DEFAULT_PROTEIN_LENGTH_BINNING       = 5
DEFAULT_GENE_CODE_ID                 = 1
# BUSCO/compleasm auto-discovery + auto-run (Python extension, not in Java)
DEFAULT_BUSCO_LINEAGE                = 'eukaryota_odb10'
DEFAULT_BUSCO_THREADS                = 4

# Regex for the bundled domains file (interproscan tabular)
_DEFAULT_DOMAINS_RE = re.compile(
    r'(?P<id>\w+)\t.*\tHMMPfam\t(?P<domain>\w+)\t.*',
    re.IGNORECASE,
)
_PFAM_IN_LINE_RE = re.compile(r'pf\d+', re.IGNORECASE)
_PROTEIN_LENGTH_RECORD_RE = re.compile(r'(\d+)\W+([0-9.\-+eE]+)')


# ---------------------------------------------------------------------------
# BUSCO / compleasm summary
#   Java port surface:  org.mycocosm.abaqs.main.BuscoData
#   Extension:          compleasm support (Java has none; see docstring header)
# ---------------------------------------------------------------------------

# Pattern for short_summary.txt file format:  "  752     Complete BUSCOs (C)"
_BUSCO_FILE_LINE_RE = re.compile(r'\s*(\d+)\s+(.+)')
# Pattern for inline string format: "C:99.3%[S:98.9%,D:0.4%],F:0.3%,M:0.4%,n:758"
_BUSCO_STRING_CAT_RE = re.compile(r'([CSDFMNn]):([0-9.]+)(%?)', re.IGNORECASE)
_BUSCO_COMMENT_RE = re.compile(r'#.*')
# compleasm summary.txt lines, e.g.  "S:97.36%, 738"
_COMPLEASM_PCT_RE = re.compile(r'\s*([SDFIM]):\s*([0-9.]+)%,\s*(\d+)')
_COMPLEASM_N_RE   = re.compile(r'\s*N:\s*(\d+)')


class BuscoData:
    """Holds BUSCO assessment proportions (all values in range 0.0–1.0)."""

    def __init__(self, complete: float, single_copy: float, duplicated: float,
                 fragmented: float, missing: float, total_groups: int):
        self.complete     = complete      # C / total
        self.single_copy  = single_copy   # S / total
        self.duplicated   = duplicated    # D / total
        self.fragmented   = fragmented    # F / total
        self.missing      = missing       # M / total
        self.total_groups = total_groups

    @staticmethod
    def of(string: str) -> 'BuscoData':
        """Parse a BUSCO summary string like
        ``C:99.3%[S:98.9%,D:0.4%],F:0.3%,M:0.4%,n:758``.

        Mirrors BuscoData.of(String) but corrects Java's percentage-only bug —
        values without a trailing ``%`` are divided by the parsed ``n`` (B1).
        """
        total_count = 0
        for m in _BUSCO_STRING_CAT_RE.finditer(string):
            if m.group(1).lower() == 'n':
                total_count = int(m.group(2))

        def _load(val_str: str, pct: str) -> float:
            if pct:
                return float(val_str) / 100.0
            if total_count > 0:
                return float(val_str) / total_count
            raise ValueError(f"Illegal BUSCO input (missing total): {string!r}")

        complete = single_copy = duplicated = fragmented = missing = 0.0
        for m in _BUSCO_STRING_CAT_RE.finditer(string):
            key = m.group(1).upper()
            val = m.group(2)
            pct = m.group(3)
            if key == 'C':
                complete    = _load(val, pct)
            elif key == 'S':
                single_copy = _load(val, pct)
            elif key == 'D':
                duplicated  = _load(val, pct)
            elif key == 'F':
                fragmented  = _load(val, pct)
            elif key == 'M':
                missing     = _load(val, pct)

        return BuscoData(complete, single_copy, duplicated,
                         fragmented, missing, total_count)

    @staticmethod
    def of_file(path: str) -> 'BuscoData':
        """Parse a BUSCO ``short_summary.txt`` file. Mirrors BuscoData.ofFile."""
        complete = single_copy = duplicated = fragmented = missing = 0
        total = 0
        with open_gzipped_or_plain(path) as fh:
            for line in fh:
                if _BUSCO_COMMENT_RE.match(line.strip()):
                    continue
                m = _BUSCO_FILE_LINE_RE.match(line)
                if not m:
                    continue
                value    = int(m.group(1))
                category = m.group(2).lower()
                if ' (c)' in category:
                    complete    = value
                elif ' (s)' in category:
                    single_copy = value
                elif ' (d)' in category:
                    duplicated  = value
                elif ' (f)' in category:
                    fragmented  = value
                elif ' (m)' in category:
                    missing     = value
                elif 'total busco ' in category:
                    total       = value

        if total == 0:
            raise ValueError(f"Could not find total BUSCO groups in file: {path!r}")

        return BuscoData(
            complete    / total,
            single_copy / total,
            duplicated  / total,
            fragmented  / total,
            missing     / total,
            total,
        )

    @staticmethod
    def of_compleasm_file(path: str) -> 'BuscoData':
        """Parse a compleasm ``summary.txt`` file.

        compleasm reports S (single-copy), D (duplicated), F (fragmented),
        I (internal-stop / incomplete), M (missing), and N (total). BUSCO's
        C = S + D; ABAQS needs C and D. compleasm's I has no BUSCO equivalent
        and is folded into ``missing`` here.

        Java has no compleasm support — this is a pure Python extension.
        """
        with open_gzipped_or_plain(path) as fh:
            return _parse_compleasm_summary_text(fh.read())


_COMPLEASM_DIR_SUFFIX = '_compleasm'
_BUSCO_DIR_SUFFIX     = '_busco'


def _pick_lineage_dir(base_dir: Path, base_name: str, suffix: str,
                      lineage_hint: Optional[str]) -> Optional[Path]:
    """Pick a `<base_name>_<lineage><suffix>` dir from *base_dir*.

    If *lineage_hint* is given and a matching dir exists, prefer it; otherwise
    return the most recently modified dir matching the glob, or None.
    """
    candidates = sorted(base_dir.glob(f"{base_name}_*{suffix}"))
    if not candidates:
        return None
    if lineage_hint:
        preferred = base_dir / f"{base_name}_{lineage_hint}{suffix}"
        if preferred in candidates:
            return preferred
    return max(candidates, key=lambda d: d.stat().st_mtime)


def _discover_busco_data(scaffolds_fasta: str,
                         lineage_hint: Optional[str]) -> Optional['BuscoData']:
    """Look for pre-existing compleasm or BUSCO output beside the FASTA.

    Convention (EGAP style): ``<base_dir>/<base_name>_<lineage>_compleasm/`` and
    ``<base_dir>/<base_name>_<lineage>_busco/``. compleasm output is preferred.
    """
    p = Path(scaffolds_fasta)
    base_dir = p.parent
    base_name = p.stem  # strips final .gz only if .gz is the last suffix; ok.

    # 1) compleasm
    chosen = _pick_lineage_dir(base_dir, base_name,
                               _COMPLEASM_DIR_SUFFIX, lineage_hint)
    if chosen is not None:
        summary = chosen / 'summary.txt'
        if summary.exists():
            log.info("Discovered compleasm summary at '%s'", summary)
            try:
                return BuscoData.of_compleasm_file(str(summary))
            except Exception as e:
                log.warning("Failed to parse compleasm summary '%s': %s", summary, e)

    # 2) BUSCO (short_summary*.txt within EGAP-style dir)
    chosen = _pick_lineage_dir(base_dir, base_name,
                               _BUSCO_DIR_SUFFIX, lineage_hint)
    if chosen is not None:
        summaries = sorted(chosen.glob('short_summary*.txt'))
        if summaries:
            target = summaries[0]
            log.info("Discovered BUSCO summary at '%s'", target)
            try:
                return BuscoData.of_file(str(target))
            except Exception as e:
                log.warning("Failed to parse BUSCO summary '%s': %s", target, e)

    return None


def _run_compleasm(scaffolds_fasta: str, output_dir: Path,
                   lineage: str, threads: int) -> bool:
    """Invoke `compleasm run` against *scaffolds_fasta*. Returns True on success."""
    import subprocess
    output_dir.mkdir(parents=True, exist_ok=True)
    cmd = ['compleasm', 'run',
           '-a', scaffolds_fasta,
           '-o', str(output_dir),
           '-l', lineage,
           '-t', str(threads)]
    log.info("Running: %s", ' '.join(cmd))
    try:
        result = subprocess.run(cmd, capture_output=True, text=True)
    except FileNotFoundError:
        log.warning("compleasm executable not found on PATH")
        return False
    except Exception as e:
        log.warning("Failed to invoke compleasm: %s", e)
        return False
    if result.returncode != 0:
        log.warning("compleasm failed (rc=%d). stderr (truncated): %s",
                    result.returncode,
                    (result.stderr or '').strip()[:500])
        return False
    return (output_dir / 'summary.txt').exists()


def _run_busco(scaffolds_fasta: str, parent_dir: Path, run_name: str,
               lineage: str, threads: int) -> bool:
    """Invoke `busco -i …` against *scaffolds_fasta*. Returns True on success."""
    import subprocess
    parent_dir.mkdir(parents=True, exist_ok=True)
    cmd = ['busco',
           '-i', scaffolds_fasta,
           '-o', run_name,
           '--out_path', str(parent_dir),
           '-l', lineage,
           '-m', 'genome',
           '-c', str(threads)]
    log.info("Running: %s", ' '.join(cmd))
    try:
        result = subprocess.run(cmd, capture_output=True, text=True)
    except FileNotFoundError:
        log.warning("busco executable not found on PATH")
        return False
    except Exception as e:
        log.warning("Failed to invoke busco: %s", e)
        return False
    if result.returncode != 0:
        log.warning("busco failed (rc=%d). stderr (truncated): %s",
                    result.returncode,
                    (result.stderr or '').strip()[:500])
        return False
    # Confirm a short_summary*.txt landed in <parent>/<run_name>/
    out_dir = parent_dir / run_name
    return bool(list(out_dir.glob('short_summary*.txt')))


def _load_or_run_busco_data(scaffolds_fasta: str, lineage: str,
                            threads: int) -> Optional['BuscoData']:
    """Discover existing compleasm/BUSCO output; else run compleasm, then BUSCO.

    Returns ``None`` if no usable result could be obtained (caller falls back
    to ``complete=1.0, duplicated=0.0`` — same as 'no BUSCO data').
    """
    found = _discover_busco_data(scaffolds_fasta, lineage)
    if found is not None:
        return found

    p = Path(scaffolds_fasta)
    base_dir, base_name = p.parent, p.stem

    # Try compleasm first
    compleasm_out = base_dir / f"{base_name}_{lineage}{_COMPLEASM_DIR_SUFFIX}"
    log.info("No existing compleasm/BUSCO output found; attempting compleasm "
             "into '%s' (lineage %s)", compleasm_out, lineage)
    if _run_compleasm(scaffolds_fasta, compleasm_out, lineage, threads):
        summary = compleasm_out / 'summary.txt'
        try:
            return BuscoData.of_compleasm_file(str(summary))
        except Exception as e:
            log.warning("compleasm ran but summary parsing failed: %s", e)

    # Fall back to BUSCO
    busco_name = f"{base_name}_{lineage}{_BUSCO_DIR_SUFFIX}"
    log.info("Falling back to BUSCO into '%s/%s'", base_dir, busco_name)
    if _run_busco(scaffolds_fasta, base_dir, busco_name, lineage, threads):
        out_dir = base_dir / busco_name
        summaries = sorted(out_dir.glob('short_summary*.txt'))
        if summaries:
            try:
                return BuscoData.of_file(str(summaries[0]))
            except Exception as e:
                log.warning("BUSCO ran but summary parsing failed: %s", e)

    log.warning("No BUSCO/compleasm data could be obtained; "
                "scoring will assume complete=1.0, duplicated=0.0")
    return None


def _parse_compleasm_summary_text(text: str) -> BuscoData:
    s = d = f = i = m = 0
    n = 0
    for line in text.splitlines():
        cm = _COMPLEASM_PCT_RE.match(line)
        if cm:
            key   = cm.group(1)
            count = int(cm.group(3))
            if   key == 'S': s = count
            elif key == 'D': d = count
            elif key == 'F': f = count
            elif key == 'I': i = count
            elif key == 'M': m = count
            continue
        nm = _COMPLEASM_N_RE.match(line)
        if nm:
            n = int(nm.group(1))
    if n == 0:
        raise ValueError("compleasm summary missing 'N:' total line")
    return BuscoData(
        complete     = (s + d) / n,
        single_copy  = s       / n,
        duplicated   = d       / n,
        fragmented   = f       / n,
        missing      = (m + i) / n,
        total_groups = n,
    )


# ---------------------------------------------------------------------------
# Data containers
# ---------------------------------------------------------------------------

@dataclass
class ProteinRecord:
    """Holds a translated / loaded protein sequence with its IDs."""
    protein_id: str
    extra:      str          # original FASTA description or mRNA ID
    sequence:   str          # amino-acid string


@dataclass
class GeneRecord:
    """Per-mRNA record holding all data needed for quality metrics.

    Mirrors org.mycocosm.abaqs.main.GeneRecord.
    """
    mrna:            Gff3Record
    mrna_count:      int               # isoform count for this gene
    cds_start:       int
    cds_end:         int
    protein:         Optional[ProteinRecord]
    domains:         Optional[Set[str]]   # normalised Pfam IDs (upper-case)
    portion_cds_masked: float
    detected_transposable_element: bool


@dataclass
class ABAQSResult:
    """All computed factors and the final score."""
    abaqs_score:               float
    protein_length_dist_factor: float
    incomplete_genes_factor:    float
    te_factor:                  float
    isoforms_factor:            float
    busco_complete_factor:      float
    busco_duplicated_factor:    float
    total_records:              int
    total_genes:                int
    total_scaffolds:            int
    total_proteins:             int
    proteins_with_domains:      int
    unique_domains:             int


# ---------------------------------------------------------------------------
# ID mapper helpers
# ---------------------------------------------------------------------------

# Templates may reference numbered groups ({0},{1},…) or named groups ({foo}).
_TEMPLATE_TOKEN_RE = re.compile(r'\{(\w+)\}')
# Java PatternTransformer flag-suffix syntax: "pattern<<CI,DA->template".
_FLAG_SPLIT_RE = re.compile(r'\s*<<\s*')
# Java GFF3 mapper definition splitter (3 colon-separated fields).
_GFF3_MAPPER_PARTS_RE = re.compile(r'\s*:\s*')

_PATTERN_FLAG_MAP = {
    'CASE_INSENSITIVE': re.IGNORECASE, 'CI': re.IGNORECASE,
    'DOTALL':           re.DOTALL,     'DA': re.DOTALL,
    'MULTILINE':        re.MULTILINE,  'ML': re.MULTILINE,
    'COMMENTS':         re.VERBOSE,    'CO': re.VERBOSE,
    'UNICODE_CASE':     re.UNICODE,    'UC': re.UNICODE,
    'UNICODE_CHARACTER_CLASS': re.UNICODE, 'UN': re.UNICODE,
    'UNIX_LINES':       0,             'UL': 0,   # no Python equivalent
    'LITERAL':          0,             'LI': 0,   # callers escape manually
    'CANON_EQ':         0,             'CE': 0,   # no Python equivalent
}


def _compile_pattern_transformer(definition: str):
    """Compile a Java PatternTransformer definition string.

    Format:   ``regex (<< FLAGS)? -> template``
    Examples: ``.+proteinId\\s*=\\s*(\\d+).*->{1}``
              ``abc(?P<n>\\d+)<<CI->prefix-{n}``

    Returns a callable ``str -> Optional[str]``. ``{0}`` is the whole match.
    Mirrors org.mycocosm.framework.text.PatternTransformer.
    """
    if '->' not in definition:
        # Java behavior: pattern only with no template → never matches usefully.
        return lambda _s: None

    pat_part, tmpl_part = definition.split('->', 1)
    pat_part  = pat_part.strip()
    tmpl_part = tmpl_part.strip()

    flags = 0
    flag_split = _FLAG_SPLIT_RE.split(pat_part, maxsplit=1)
    if len(flag_split) == 2:
        pat_part = flag_split[0].strip()
        for tok in re.split(r'[\s;,]+', flag_split[1].strip().upper()):
            if tok in _PATTERN_FLAG_MAP:
                flags |= _PATTERN_FLAG_MAP[tok]

    # Translate Java named-group syntax "(?<name>" to Python's "(?P<name>"
    # so users can copy Java patterns verbatim. Leaves Python-native (?P<…)
    # untouched (P would already be present).
    pat_part = re.sub(r'\(\?<([A-Za-z_][\w]*)>', r'(?P<\1>', pat_part)

    compiled = re.compile(pat_part, flags)

    def _resolve(m: re.Match) -> str:
        def _sub(tok: re.Match) -> str:
            key = tok.group(1)
            try:
                return m.group(int(key) if key.isdigit() else key) or ''
            except (IndexError, re.error):
                return ''
        return _TEMPLATE_TOKEN_RE.sub(_sub, tmpl_part)

    def _apply(s: str) -> Optional[str]:
        if s is None:
            return None
        m = compiled.fullmatch(s)
        return _resolve(m) if m else None

    return _apply


def _make_gff3_protein_id_mapper(definition: str = DEFAULT_GFF3_PROTEIN_ID_ATTR):
    """Build a GFF3 record → protein-ID mapper.

    Accepts both:
      - A bare attribute name (Python convenience): ``"proteinId"``
      - The Java MapperDefinitionsFactory format: ``"column:attribute:pattern->template"``
        where column ∈ {seqid, source, attributes}.

    Mirrors org.mycocosm.abaqs.main.MapperDefinitionsFactory.Gff3RecordIdMapper.
    """
    if ':' in definition:
        parts = _GFF3_MAPPER_PARTS_RE.split(definition, maxsplit=2)
        column     = parts[0] if len(parts) >= 1 else 'attributes'
        attr_name  = parts[1] if len(parts) >= 2 else 'proteinId'
        pat_def    = parts[2] if len(parts) >= 3 else None
    else:
        column, attr_name, pat_def = 'attributes', definition, None

    if column not in ('seqid', 'source', 'attributes'):
        raise ValueError(f"Invalid GFF3 column in mapper '{definition}': {column!r}")

    transformer = _compile_pattern_transformer(pat_def) if pat_def else None

    def _mapper(rec: Gff3Record) -> Optional[str]:
        if rec is None:
            return None
        if column == 'seqid':
            value = rec.seqid
        elif column == 'source':
            value = rec.source
        else:
            value = rec.attributes.get(attr_name)
        if value is None or value == '':
            return value
        return transformer(value) if transformer else value

    return _mapper


def _make_fasta_protein_id_mapper(definition: str):
    """Build a FASTA-header → protein-ID mapper.

    Uses Java PatternTransformer syntax with support for numbered ({0},{1},…)
    AND named (?P<n>…) → {n} groups, plus optional ``<<CI,DA`` flag suffix.
    Input passed to the regex is ``fasta_id + ' ' + extra`` (Java behavior).
    Default: ``.+proteinId\\s*=\\s*(\\d+).*->{1}``
    """
    transformer = _compile_pattern_transformer(definition)

    def _mapper(fasta_id: str, extra: str) -> Optional[str]:
        combined = (fasta_id + ' ' + (extra or '')).strip()
        return transformer(combined)

    return _mapper


# ---------------------------------------------------------------------------
# File loaders
# ---------------------------------------------------------------------------

def _load_scaffolds_fasta(path: str) -> Dict[str, str]:
    log.info("Loading scaffolds FASTA from '%s'", path)
    scaffolds = load_fasta(path)
    log.info("Loaded %d scaffolds", len(scaffolds))
    return scaffolds


def _load_proteins_fasta(path: str, id_mapper) -> Dict[str, ProteinRecord]:
    """Load a protein FASTA file; return dict protein_id → ProteinRecord.

    *id_mapper(fasta_id, extra)* → protein_id string or None.
    Mirrors ABAQS.loadProteinsFastaWithIdMapper().
    """
    log.info("Loading proteins FASTA from '%s'", path)
    proteins: Dict[str, ProteinRecord] = {}
    current_id: Optional[str] = None
    current_extra: str = ''
    current_seq: list = []

    def _flush():
        if current_id is None:
            return
        seq = ''.join(current_seq)
        pid = id_mapper(current_id, current_extra)
        if pid:
            proteins[pid] = ProteinRecord(pid, current_extra, seq)

    import gzip as _gzip
    opener = _gzip.open if path.endswith('.gz') else open
    with opener(path, 'rt') as fh:
        for line in fh:
            line = line.rstrip('\n\r')
            if line.startswith('>'):
                _flush()
                parts = line[1:].split(None, 1)
                current_id    = parts[0] if parts else ''
                current_extra = parts[1] if len(parts) > 1 else ''
                current_seq   = []
            else:
                current_seq.append(line)
    _flush()
    log.info("Loaded %d proteins", len(proteins))
    return proteins


def _translate_proteins(gff3: GFF3Data, id_mapper, gene_code: GeneCode,
                        verbose: bool = False) -> Dict[str, ProteinRecord]:
    """Translate CDS from GFF3 using *gene_code*."""
    log.info("Translating proteins using gene code '%s' (%d)",
             gene_code.names[0], gene_code.id)
    proteins: Dict[str, ProteinRecord] = {}
    for mrna in gff3.get_records_by_predicate(lambda r: r.type == Gff3Type.mRNA):
        scaffold_seq = gff3.scaffolds.get(mrna.seqid)
        if scaffold_seq is None:
            continue
        seq, overhang = mrna.get_translated_aminoacid_sequence(scaffold_seq, gene_code)
        pid = id_mapper(mrna)
        if pid is None:
            continue
        if verbose and overhang:
            log.warning("GFF record '%s' has overhang after translation: '%s'",
                        mrna.id, overhang)
        proteins[pid] = ProteinRecord(pid, mrna.id, seq)
    return proteins


def _load_domains(path: Optional[str],
                  domains_re: re.Pattern,
                  verbose: bool = False) -> Dict[str, Set[str]]:
    """Load Pfam domain assignments; return dict protein_id → set(pfam_ids)."""
    domains: Dict[str, Set[str]] = {}
    if path is None:
        log.warning("No domains file specified")
        return domains
    log.info("Loading domains from '%s'", path)
    skipped = added = 0
    with open_gzipped_or_plain(path) as fh:
        for line in fh:
            line = line.rstrip('\n\r')
            m = domains_re.fullmatch(line)
            if m:
                pid   = m.group('id')
                pfam  = m.group('domain').upper()
                if pid not in domains:
                    domains[pid] = set()
                domains[pid].add(pfam)
                added += 1
                if verbose:
                    log.info("Added domain '%s' to protein '%s'", pfam, pid)
            else:
                skipped += 1
                if verbose:
                    log.warning("Skipped domain line: '%s'", line)
    log.info("Added %d domains to %d proteins; skipped %d lines",
             added, len(domains), skipped)
    return domains


def _load_te_domains(path: Optional[str], resource_name: str) -> Set[str]:
    """Load TE Pfam domain list (file or bundled resource); return normalised set."""
    if path is not None:
        log.info("Loading TE domains from file '%s'", path)
        opener = gzip.open if path.endswith('.gz') else open
        with opener(path, 'rt') as fh:
            text = fh.read()
    else:
        log.info("Using bundled %s", resource_name)
        data_dir = Path(__file__).parent / 'data'
        fpath = data_dir / resource_name
        text = fpath.read_text(encoding='utf-8', errors='replace')

    domains: Set[str] = set()
    for line in text.splitlines():
        line = line.strip()
        if not line or line.startswith('#'):
            continue
        for m in _PFAM_IN_LINE_RE.finditer(line):
            domains.add(m.group(0).upper())
    log.info("Loaded %d %s domains", len(domains), resource_name)
    return domains


def _load_reference_protein_length_distribution(path: Optional[str]) -> Dict[int, float]:
    """Load and normalise the reference protein-length distribution.

    Keys are the integers in the first column (bin index or protein length).
    Values are normalised so they sum to 1.0.
    Mirrors ABAQS.loadReferenceProteinLengthDistribution().
    """
    if path is not None:
        log.info("Loading reference distribution from '%s'", path)
        opener = gzip.open if path.endswith('.gz') else open
        with opener(path, 'rt') as fh:
            text = fh.read()
    else:
        log.info("Using bundled reference protein-length distribution")
        data_dir = Path(__file__).parent / 'data'
        with gzip.open(str(data_dir / 'reference-proteins-length-distribution.tsv.gz'), 'rt') as fh:
            text = fh.read()

    raw: Dict[int, float] = {}
    for line in text.splitlines():
        m = _PROTEIN_LENGTH_RECORD_RE.fullmatch(line.strip())
        if m:
            raw[int(m.group(1))] = float(m.group(2))

    total = sum(raw.values())
    if total == 0.0:
        return raw
    return {k: v / total for k, v in raw.items()}


# ---------------------------------------------------------------------------
# Factor computations
# ---------------------------------------------------------------------------

def _compute_cds_masked_factor(scaffold_seq: str, mrna: Gff3Record,
                                masker: str) -> float:
    """Fraction of CDS bases that are soft-masked.

    Mirrors ABAQS.computeScaffoldMaskedFactor().
    """
    total = masked = 0
    for seg_start, seg_end in mrna.get_segments_by_type(Gff3Type.CDS):
        for pos in range(seg_start, seg_end + 1):
            ch = scaffold_seq[pos - 1]
            total += 1
            if masker == 'TO_LOWER_CASE' and ch.islower():
                masked += 1
            elif masker == 'TO_N' and ch.upper() == 'N':
                masked += 1
            elif masker == 'TO_DOT' and ch == '.':
                masked += 1
    return masked / total if total else 0.0


def _is_transposable_element(
    domains: Optional[Set[str]],
    portion_cds_masked: float,
    te_domains: Set[str],
    suspected_te_domains: Set[str],
    no_domain_cutoff: float,
    suspected_cutoff: float,
) -> bool:
    """Classify a gene model as a transposable element.

    Mirrors ABAQS.isTransposableElement().
    Rules (in order):
      a. Any domain in known TE list → TE
      c. Suspected TE domain AND (cutoff is NaN OR masked > cutoff) → TE
    Rule (b) from the Java source — "no domains AND CDS masking > cutoff" —
    is commented out in Java and not implemented here. The
    ``no_domain_cds_masked_cutoff`` parameter is accepted but unused.
    """
    if domains:
        # Rule a
        if any(d in te_domains for d in domains):
            return True
        # Rule c (suspected)
        if any(d in suspected_te_domains for d in domains):
            if math.isnan(suspected_cutoff) or portion_cds_masked > suspected_cutoff:
                return True
    return False


def _create_gene_records(
    gff3: GFF3Data,
    masker: str,
    proteins: Dict[str, ProteinRecord],
    domains: Dict[str, Set[str]],
    gff3_id_mapper,
    te_domains: Set[str],
    suspected_te_domains: Set[str],
    no_domain_cutoff: float,
    suspected_cutoff: float,
) -> Dict[str, GeneRecord]:
    """Build per-mRNA GeneRecord objects.

    Mirrors ABAQS.createGeneRecords().
    """
    records: Dict[str, GeneRecord] = {}

    for gene_rec in gff3.records:
        if gene_rec.category != Gff3RecordCategory.regular:
            continue
        mrnas = gene_rec.get_all_by_predicate(lambda r: r.type == Gff3Type.mRNA)

        # Gene-level CDS span (across all isoforms, matching Java behaviour)
        all_cds = gene_rec.get_all_by_predicate(lambda r: r.type == Gff3Type.CDS)
        cds_start = min((c.start for c in all_cds), default=0)
        cds_end   = max((c.end   for c in all_cds), default=0)

        for mrna in mrnas:
            pid = gff3_id_mapper(mrna)
            scaffold_seq = gff3.scaffolds.get(mrna.seqid)
            masked = (
                _compute_cds_masked_factor(scaffold_seq, mrna, masker)
                if scaffold_seq else 0.0
            )
            dom = domains.get(pid) if pid else None
            is_te = _is_transposable_element(
                dom, masked, te_domains, suspected_te_domains,
                no_domain_cutoff, suspected_cutoff,
            )
            prot = proteins.get(pid) if pid else None
            records[mrna.id] = GeneRecord(
                mrna            = mrna,
                mrna_count      = len(mrnas),
                cds_start       = cds_start,
                cds_end         = cds_end,
                protein         = prot,
                domains         = dom,
                portion_cds_masked = masked,
                detected_transposable_element = is_te,
            )

    return records


def _compute_incomplete_genes_factor(gene_records: Dict[str, GeneRecord]) -> float:
    """Fraction of genes that start with M and end with *.

    Returns 1 - incomplete_fraction.
    Mirrors ABAQS.computeIncompleteGenesFactor().
    """
    incomplete = 0
    for gr in gene_records.values():
        if gr.protein and gr.protein.sequence:
            seq = gr.protein.sequence
            if seq[0] != 'M' or seq[-1] != '*':
                incomplete += 1
    log.info("Found %d incomplete genes", incomplete)
    return 1.0 - incomplete / len(gene_records) if gene_records else 1.0


def _compute_te_factor(gene_records: Dict[str, GeneRecord]) -> float:
    """TE count / total gene count.

    Mirrors ABAQS.computeTransposableElementsFactor().
    """
    te_count = sum(1 for gr in gene_records.values()
                   if gr.detected_transposable_element)
    log.info("Found %d TE/suspected TE", te_count)
    return te_count / len(gene_records) if gene_records else 0.0


def _compute_isoforms_factor(
    gff3: GFF3Data,
    gene_records: Dict[str, GeneRecord],
    isoforms_min_overlap: float,
) -> float:
    """Isoform factor: (gff3_isoform_genes + discovered_overlapping) / total_mRNAs.

    Mirrors ABAQS.computeIsoformsFactor().
    """
    track = FeatureTrack.create_from_gff3(gff3, 'gff3Data', model=True)
    self_mapping = FeatureTrack.map_two_tracks_by_position_overlap(
        track, track, FeaturePair.compare_by_overlap_desc
    )

    # GFF3-defined isoforms: genes with > 1 mRNA child
    gff3_isoform_mrna_ids: Set[str] = set()
    isoforms_count = 0
    for gene in gff3.get_records_by_predicate(lambda r: r.type == Gff3Type.gene):
        mrnas = gene.get_all_by_predicate(lambda r: r.type == Gff3Type.mRNA)
        if len(mrnas) > 1:
            isoforms_count += 1
            for m in mrnas:
                gff3_isoform_mrna_ids.add(m.id)

    log.info("Found %d GFF3 isoform genes", isoforms_count)

    # Discovered overlapping pairs not already flagged as GFF3 isoforms
    overlapped_count = 0
    for pair in self_mapping.mapped_pairs:
        if pair.overlap > isoforms_min_overlap:
            id_a = pair.a.name
            id_b = pair.b.name
            if id_a not in gff3_isoform_mrna_ids or id_b not in gff3_isoform_mrna_ids:
                overlapped_count += 1
                gff3_isoform_mrna_ids.add(id_a)
                gff3_isoform_mrna_ids.add(id_b)

    log.info("Found %d genes overlapped by > %.0f%%",
             overlapped_count, isoforms_min_overlap * 100)

    total = len(gene_records)
    return (isoforms_count + overlapped_count) / total if total else 0.0


def _compute_organism_protein_length_distribution(
    proteins: list,   # list of ProteinRecord (non-TE)
    bin_size: int,
) -> Dict[int, float]:
    """Bin organism proteins by length; normalise counts.

    Key = bin index (0, 1, 2, …).
    Mirrors ABAQS.loadOrganismProteinLengthDistribution().
    """
    if not proteins:
        return {}
    total = len(proteins)
    binned = bin_collection_fixed_bin_size(
        proteins,
        to_bin_value=lambda p: float(len(p.sequence)),
        bin_step=float(bin_size),
        min_value=0.0,
        max_value=None,
    )
    result: Dict[int, float] = {}
    for idx, (_, items) in enumerate(binned):
        result[idx] = len(items) / total
    log.info("Binned %d proteins into %d bins", total, len(result))
    return result


def _compute_protein_length_distribution_factor(
    organism_dist: Dict[int, float],
    reference_dist: Dict[int, float],
) -> float:
    """Compare organism vs. reference protein-length distributions.

    Formula:
        diff  = (length + 0.5) * (ref_value - org_value)   for each shared key
        base  = (length + 0.5) * ref_value
        factor = clamp(1.0 - sum(diff) / sum(base), 0.0, 1.0)

    Mirrors ABAQS.computeProteinLengthDistributionFactor().
    """
    all_keys = set(organism_dist) | set(reference_dist)
    sum_diff = 0.0
    sum_base = 0.0
    for k in all_keys:
        if k <= 0:
            continue
        org_val = organism_dist.get(k, 0.0)
        ref_val = reference_dist.get(k, 0.0)
        weight   = k + 0.5
        sum_diff += weight * (ref_val - org_val)
        sum_base += weight * ref_val

    if sum_base == 0.0:
        return 1.0
    raw = 1.0 - sum_diff / sum_base
    return max(0.0, min(1.0, raw))


# ---------------------------------------------------------------------------
# Main public API
# ---------------------------------------------------------------------------

def compute_abaqs(
    input_gff3:                str,
    input_scaffolds_fasta:     Optional[str]  = None,
    input_proteins_fasta:      Optional[str]  = None,
    input_domains:             Optional[str]  = None,
    busco_string:              Optional[str]  = None,
    busco_file:                Optional[str]  = None,
    te_domains_file:           Optional[str]  = None,
    suspected_te_domains_file: Optional[str]  = None,
    reference_pld_file:        Optional[str]  = None,
    gene_code_file:            Optional[str]  = None,
    gene_code_id:              int            = DEFAULT_GENE_CODE_ID,
    gff3_protein_id_attr:      str            = DEFAULT_GFF3_PROTEIN_ID_ATTR,
    protein_fasta_id_pattern:  Optional[str]  = None,
    domains_pattern:           Optional[str]  = None,
    isoforms_min_overlap:      float          = DEFAULT_ISOFORMS_MIN_OVERLAP,
    masker_function:           str            = DEFAULT_MASKER_FUNCTION,
    no_domain_cds_masked_cutoff: float        = DEFAULT_NO_DOMAINS_CDS_MASKED_CUTOFF,
    suspected_domain_masked_cutoff: float     = DEFAULT_SUSPECTED_DOMAINS_MASKED_CUTOFF,
    protein_length_binning:    int            = DEFAULT_PROTEIN_LENGTH_BINNING,
    busco_auto_run:            bool           = True,
    busco_lineage:             str            = DEFAULT_BUSCO_LINEAGE,
    busco_threads:             int            = DEFAULT_BUSCO_THREADS,
    verbose:                   bool           = False,
) -> ABAQSResult:
    """Run the full ABAQS scoring pipeline and return an ABAQSResult.

    Mirrors ABAQS.processInput().
    """
    log.info("Starting ABAQS v1.0")

    # ---- Gene code ----
    if gene_code_file:
        codes = load_gene_codes_file(gene_code_file)
    else:
        codes = load_default_gene_codes()
    gene_code = codes[gene_code_id]
    log.info("Loaded gene code '%s' (%d)", gene_code.names[0], gene_code.id)

    # ---- GFF3 ----
    log.info("Loading GFF3 from '%s'", input_gff3)
    gff3 = parse_gff3_file(input_gff3)
    total_records = len(gff3.records)
    total_genes   = sum(1 for r in gff3.records if r.category.value == 'regular' and r.parent is None)
    log.info("Loaded %d records (%d top-level genes)", total_records, total_genes)

    # ---- Scaffolds ----
    if input_scaffolds_fasta:
        scaffolds = _load_scaffolds_fasta(input_scaffolds_fasta)
        gff3 = gff3.replace_scaffolds(scaffolds)
    elif not gff3.has_scaffolds():
        raise ValueError(
            "Must supply either --input-scaffolds-fasta or a GFF3 with embedded ##FASTA"
        )
    total_scaffolds = len(gff3.scaffolds)

    # ---- ID mappers ----
    gff3_id_mapper = _make_gff3_protein_id_mapper(gff3_protein_id_attr)
    if protein_fasta_id_pattern:
        fasta_id_mapper = _make_fasta_protein_id_mapper(protein_fasta_id_pattern)
    else:
        # Default: ".+proteinId\s*=\s*(\d+).*->{1}"
        fasta_id_mapper = _make_fasta_protein_id_mapper(
            r'.+proteinId\s*=\s*(\d+).*->{1}'
        )

    # ---- Proteins ----
    if input_proteins_fasta:
        proteins_dict = _load_proteins_fasta(input_proteins_fasta, fasta_id_mapper)
    else:
        proteins_dict = _translate_proteins(gff3, gff3_id_mapper, gene_code, verbose)

    total_proteins = len(proteins_dict)

    # ---- Domains ----
    if domains_pattern:
        dom_re = re.compile(domains_pattern, re.IGNORECASE)
    else:
        dom_re = _DEFAULT_DOMAINS_RE
    domains = _load_domains(input_domains, dom_re, verbose)

    # ---- TE domain lists ----
    te_domains          = _load_te_domains(te_domains_file,
                                           'transposable-elements-pfams.txt')
    suspected_te_domains = _load_te_domains(suspected_te_domains_file,
                                             'suspected-transposable-elements-pfams.txt')

    # Normalise domain keys to upper-case for comparison
    te_domains_upper          = {d.upper() for d in te_domains}
    suspected_te_domains_upper = {d.upper() for d in suspected_te_domains}

    # ---- GeneRecords ----
    gene_records = _create_gene_records(
        gff3, masker_function, proteins_dict, domains,
        gff3_id_mapper, te_domains_upper, suspected_te_domains_upper,
        no_domain_cds_masked_cutoff, suspected_domain_masked_cutoff,
    )
    log.info("Built %d gene records", len(gene_records))

    # ---- Isoforms factor ----
    isoforms_factor = _compute_isoforms_factor(gff3, gene_records, isoforms_min_overlap)
    log.info("Isoforms factor: %.4f", isoforms_factor)

    # ---- BUSCO / compleasm ----
    # Resolution order (Python extension over Java):
    #   1. Explicit --busco-data string         (-ib)
    #   2. Explicit --busco-data-file path      (-ibf)
    #   3. Auto-discover beside the scaffolds FASTA: compleasm summary first,
    #      then BUSCO summary, matching EGAP-style
    #      '<basename>_<lineage>_compleasm/summary.txt' or
    #      '<basename>_<lineage>_busco/short_summary*.txt'.
    #   4. Run compleasm; if that fails, run BUSCO.
    #   5. Warn and skip (downstream sets complete=1.0, duplicated=0.0).
    if busco_string:
        busco = BuscoData.of(busco_string)
        log.info("Parsed BUSCO from string: '%s'", busco_string)
    elif busco_file:
        # Allow either format: peek at extension; compleasm output is
        # conventionally named 'summary.txt', BUSCO 'short_summary*.txt'.
        bf_name = Path(busco_file).name.lower()
        if bf_name == 'summary.txt' or 'compleasm' in bf_name:
            busco = BuscoData.of_compleasm_file(busco_file)
            log.info("Loaded compleasm summary from file: '%s'", busco_file)
        else:
            busco = BuscoData.of_file(busco_file)
            log.info("Loaded BUSCO short_summary from file: '%s'", busco_file)
    elif busco_auto_run and input_scaffolds_fasta:
        busco = _load_or_run_busco_data(
            input_scaffolds_fasta, busco_lineage, busco_threads)
    else:
        busco = None
        log.info("No BUSCO/compleasm data provided and auto-run disabled")

    if busco is not None:
        busco_complete_factor    = busco.complete
        busco_duplicated_factor  = busco.duplicated
    else:
        busco_complete_factor    = 1.0
        busco_duplicated_factor  = 0.0
    log.info("BUSCO complete factor: %.4f", busco_complete_factor)
    log.info("BUSCO duplicated factor: %.4f", busco_duplicated_factor)

    # ---- Incomplete genes factor ----
    incomplete_factor = _compute_incomplete_genes_factor(gene_records)
    log.info("Incomplete genes factor: %.4f", incomplete_factor)

    # ---- TE factor ----
    te_factor = _compute_te_factor(gene_records)
    log.info("TE factor: %.4f", te_factor)

    # ---- Protein length distribution factor ----
    non_te_proteins = [
        gr.protein
        for gr in gene_records.values()
        if not gr.detected_transposable_element and gr.protein is not None
    ]
    org_pld = _compute_organism_protein_length_distribution(
        non_te_proteins, protein_length_binning
    )
    ref_pld = _load_reference_protein_length_distribution(reference_pld_file)
    pld_factor = _compute_protein_length_distribution_factor(org_pld, ref_pld)
    log.info("Protein length distribution factor: %.4f", pld_factor)

    # ---- Final ABAQS score ----
    upper      = math.sqrt(pld_factor * incomplete_factor)
    lower      = 1.0 + 0.5 * (te_factor + isoforms_factor +
                               busco_duplicated_factor + (1.0 - busco_complete_factor))
    abaqs_score = upper / lower
    log.info("ABAQS score: %.4f", abaqs_score)

    # ---- Stats for output ----
    proteins_with_domains = sum(
        1 for pid in proteins_dict if pid in domains
    )
    unique_domains: set = set()
    for pid, dom_set in domains.items():
        if pid in proteins_dict:
            unique_domains |= dom_set

    return ABAQSResult(
        abaqs_score               = abaqs_score,
        protein_length_dist_factor = pld_factor,
        incomplete_genes_factor    = incomplete_factor,
        te_factor                  = te_factor,
        isoforms_factor            = isoforms_factor,
        busco_complete_factor      = busco_complete_factor,
        busco_duplicated_factor    = busco_duplicated_factor,
        total_records              = total_records,
        total_genes                = total_genes,
        total_scaffolds            = total_scaffolds,
        total_proteins             = total_proteins,
        proteins_with_domains      = proteins_with_domains,
        unique_domains             = len(unique_domains),
    )


def print_result(result: ABAQSResult, gff3_path: str, out=None) -> None:
    """Print the ABAQS result in the same tab-separated format as the Java tool."""
    import sys
    fh = out or sys.stdout
    fh.write(f"Computing ABAQS score from input:\t'{gff3_path}'\n")
    fh.write(f"Total records:\t{result.total_records}\n")
    fh.write(f"Total genes:\t{result.total_genes}\n")
    fh.write(f"Total scaffolds:\t{result.total_scaffolds}\n")
    fh.write(f"Total proteins:\t{result.total_proteins}\n")
    fh.write(f"Total proteins with domains:\t{result.proteins_with_domains}\n")
    fh.write(f"Total unique domains:\t{result.unique_domains}\n")
    fh.write(f"Protein lengths distrbution factor:\t{result.protein_length_dist_factor:.4f}\n")
    fh.write(f"Incomplete genes factor:\t{result.incomplete_genes_factor:.4f}\n")
    fh.write(f"Transposable elements factor:\t{result.te_factor:.4f}\n")
    fh.write(f"Isoforms factor:\t{result.isoforms_factor:.4f}\n")
    fh.write(f"BUSCO duplicated factor:\t{result.busco_duplicated_factor:.4f}\n")
    fh.write(f"BUSCO complete factor:\t{result.busco_complete_factor:.4f}\n")
    fh.write(f"ABAQS score:\t{result.abaqs_score:.4f}\n")
