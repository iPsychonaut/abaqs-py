"""Biological sequence parsing and transformation utilities.

Consolidates three former modules:
  - GFF3 record/data parsing      (was ``abaqs.gff3``)
  - NCBI genetic codes/translation (was ``abaqs.ncbi``)
  - Fixed-bin-size binning         (was ``abaqs.binning``)

Mirrors:
  org.mycocosm.gff3.{Gff3Record,GFF3Data,Gff3Type,Gff3RecordCategory}
  org.mycocosm.ncbi.{GeneCode,GeneCodeHelper}
  org.mycocosm.framework.biology.{Strand,DnaHelper}
  org.mycocosm.framework.fasta.* (the bits used by ABAQS)
  org.mycocosm.framework.collections.CollectionsHelper#binCollectionWithFixedBinSize

Intentional divergences from Java:
  - A5: ``Strand.of`` returns ``Strand.unknown`` for unrecognized characters
        (logs a warning) instead of throwing IllegalArgumentException; empty
        input returns ``Strand.unknown`` instead of ``null``.
  - A5: ``Gff3Record.from_line`` wraps parse failures with the original line
        as context (matches Java's ExceptionsHelper wrapping).
  - C4: ``load_gene_codes_file(path)`` and ``parse_gene_codes_text(text)``
        are separate entry points (Java takes a single InputStream).
  - D3: FASTA reader accepts ``>`` and ``;`` header lines; Java only ``>``.
  - Binning: empty input returns ``[]`` instead of indexing an empty list.
"""

from __future__ import annotations

import gzip
import importlib.resources as _ir
import logging
import math
import re
from enum import Enum
from pathlib import Path
from typing import Callable, Dict, List, Optional, Tuple, TypeVar

log = logging.getLogger(__name__)


# =============================================================================
# Section 1 — Enumerations
# =============================================================================

class Strand(Enum):
    plus    = '+'
    minus   = '-'
    none    = '.'
    unknown = '?'

    @staticmethod
    def of(s: str) -> 'Strand':
        if not s:
            return Strand.unknown
        ch = s[0]
        mapped = {'+': Strand.plus,
                  '-': Strand.minus,
                  '.': Strand.none,
                  '?': Strand.unknown}.get(ch)
        if mapped is None:
            log.warning("Unrecognized strand character %r; coercing to '?'", ch)
            return Strand.unknown
        return mapped


class Gff3Type(Enum):
    gene                  = 'gene'
    mRNA                  = 'mRNA'
    tRNA                  = 'tRNA'
    rRNA                  = 'rRNA'
    five_prime_UTR        = 'five_prime_UTR'
    three_prime_UTR       = 'three_prime_UTR'
    exon                  = 'exon'
    CDS                   = 'CDS'
    intron                = 'intron'
    TF_binding_site       = 'TF_binding_site'
    mobile_genetic_element = 'mobile_genetic_element'
    repeat_region         = 'repeat_region'

    @staticmethod
    def value_of_or_none(s: str) -> Optional['Gff3Type']:
        try:
            return Gff3Type(s)
        except ValueError:
            return None


class Gff3RecordCategory(Enum):
    regular = 'regular'
    meta    = 'meta'
    fasta   = 'fasta'


# =============================================================================
# Section 2 — URL-encoding helpers (GFF3 spec)
# =============================================================================

_ESCAPED_RE = re.compile(r'%([a-fA-F0-9]{2})')
_MUST_ESCAPE = set('\t\n\r%;=&,\x7f') | set(chr(c) for c in range(0x20))


def _unescape(s: str) -> str:
    if not s:
        return s
    return _ESCAPED_RE.sub(lambda m: chr(int(m.group(1), 16)), s)


def _escape(s: str) -> str:
    if s is None:
        return ''
    result = []
    for ch in s:
        if ch in _MUST_ESCAPE or ord(ch) < 0x20:
            result.append(f'%{ord(ch):02X}')
        else:
            result.append(ch)
    return ''.join(result)


# =============================================================================
# Section 3 — Low-level column parsers
# =============================================================================

_COMMENT_RE      = re.compile(r'#\s+.*|#$')
_META_RE         = re.compile(r'##(\S+)\s*(.*)')
_POSITION_RE     = re.compile(r'[<>]?(\d+)')
_ATTR_PARSE_RE   = re.compile(r'(\S+)=(.*)')


def _parse_position(s: str) -> int:
    m = _POSITION_RE.fullmatch(s.strip())
    if m:
        return int(m.group(1))
    raise ValueError(f"Invalid position: {s!r}")


def _parse_score(s: str) -> float:
    return float('nan') if s in ('.', '') else float(s)


def _parse_phase(s: str) -> int:
    return -1 if s in ('.', '') else int(s)


def _parse_attributes(attr_str: str) -> Dict[str, str]:
    if not attr_str or attr_str == '.':
        return {}
    result: Dict[str, str] = {}
    for part in attr_str.split(';'):
        part = part.strip()
        if not part:
            continue
        m = _ATTR_PARSE_RE.match(part)
        if m:
            result[_unescape(m.group(1))] = _unescape(m.group(2))
        else:
            raise ValueError(f"Cannot parse GFF3 attribute: {part!r}")
    return result


# =============================================================================
# Section 4 — DNA helpers
# =============================================================================

_COMP_TABLE = str.maketrans('ATCGatcg', 'TAGCtagc')


def _reverse_complement(seq: str) -> str:
    return seq.translate(_COMP_TABLE)[::-1]


# =============================================================================
# Section 5 — NCBI genetic codes
# =============================================================================

DEFAULT_GENE_CODE      = 1  # Standard
DEFAULT_MITO_GENE_CODE = 4  # Vertebrate mitochondrial

_VALID_NUCLEOTIDES = re.compile(r'[ATCG]+')


class GeneCode:
    """A single NCBI genetic code translation table.

    Mirrors org.mycocosm.ncbi.GeneCode.
    """

    def __init__(self, id_: int, names: list,
                 ncbieaa: Dict[str, str],
                 sncbieaa: Dict[str, str],
                 base_1: str, base_2: str, base_3: str):
        self.id       = id_
        self.names    = names
        self.ncbieaa  = ncbieaa   # codon → amino-acid character
        self.sncbieaa = sncbieaa  # codon → start codon character
        self.base_1   = base_1
        self.base_2   = base_2
        self.base_3   = base_3

    def translate_to_amino_acid(self, nucleotide_seq: str) -> Tuple[str, Optional[str]]:
        return translate_sequence(nucleotide_seq, self)


def translate_sequence(sequence: str, code: GeneCode) -> Tuple[str, Optional[str]]:
    """Translate *sequence* using *code*.

    Returns ``(output_protein, overhang)`` where ``overhang`` is ``None`` or the
    1-2 trailing bases that don't form a complete codon. Non-ATCG triplets
    become ``'x'``. Mirrors ``GeneCode.translate()``.
    """
    result = []
    idx = 0
    table = code.ncbieaa
    while idx <= len(sequence) - 3:
        triplet = sequence[idx:idx + 3].upper()
        if _VALID_NUCLEOTIDES.fullmatch(triplet):
            aa = table.get(triplet)
            if aa is None:
                raise ValueError(
                    f"Missing triplet '{triplet}' in translation table "
                    f"'{code.names[0]}' at position {idx}"
                )
            result.append(aa)
        else:
            result.append('x')
        idx += 3

    overhang = sequence[idx:] if idx < len(sequence) else None
    return ''.join(result), overhang


# ---- gc.prt parser --------------------------------------------------------

_GC_SECTION_RE = re.compile(r'\{([^}]+)\}', re.DOTALL | re.MULTILINE)
_GC_SPLIT_LINE_RE = re.compile(r'[\n\r\t{,]+')
_GC_BASE_START_RE = re.compile(r'^[\s\-]+')
_GC_KEYVAL_RE     = re.compile(r'(\w+)\s+(.+)')


def _trim_quotes(s: str, *chars) -> str:
    for ch in chars:
        s = s.strip(ch)
    return s


def _parse_gc_section(section_text: str) -> Optional[GeneCode]:
    names: list  = []
    id_    = -1
    ncbieaa_str  = None
    sncbieaa_str = None
    base_1 = base_2 = base_3 = None

    for line in _GC_SPLIT_LINE_RE.split(section_text):
        line = line.strip()
        if not line:
            continue
        line = _GC_BASE_START_RE.sub('', line)
        m = _GC_KEYVAL_RE.match(line)
        if not m:
            continue
        key   = m.group(1).lower()
        value = _trim_quotes(m.group(2).strip(), ' ', '"')
        if key == 'name':
            names.append(value)
        elif key == 'id':
            id_ = int(value)
        elif key == 'ncbieaa':
            ncbieaa_str = value
        elif key == 'sncbieaa':
            sncbieaa_str = value
        elif key == 'base1':
            base_1 = value
        elif key == 'base2':
            base_2 = value
        elif key == 'base3':
            base_3 = value

    if id_ < 0 or base_1 is None:
        return None

    ncbieaa_map:  Dict[str, str] = {}
    sncbieaa_map: Dict[str, str] = {}
    for i in range(len(base_1)):
        codon = base_1[i] + base_2[i] + base_3[i]
        ncbieaa_map[codon]  = ncbieaa_str[i]
        sncbieaa_map[codon] = sncbieaa_str[i]

    return GeneCode(id_, names, ncbieaa_map, sncbieaa_map, base_1, base_2, base_3)


def parse_gene_codes_text(text: str) -> Dict[int, GeneCode]:
    """Parse NCBI gc.prt content; return mapping of code-id → GeneCode.

    Mirrors GeneCodeHelper.loadAndParseGeneCodes() but takes already-decoded text.
    """
    codes: Dict[int, GeneCode] = {}
    for m in _GC_SECTION_RE.finditer(text):
        code = _parse_gc_section(m.group(1))
        if code is not None:
            codes[code.id] = code
    return codes


def load_gene_codes_file(path: str) -> Dict[int, GeneCode]:
    """Load and parse an NCBI gc.prt file (plain or gzipped)."""
    p = Path(path)
    opener = gzip.open if p.suffix == '.gz' else open
    with opener(str(p), 'rt') as fh:
        return parse_gene_codes_text(fh.read())


def load_gene_codes(path_or_text: str) -> Dict[int, GeneCode]:
    """Compatibility wrapper: dispatch on whether the argument looks like a path.

    Prefer ``load_gene_codes_file`` or ``parse_gene_codes_text`` in new code.
    """
    if len(path_or_text) < 512 and Path(path_or_text).exists():
        return load_gene_codes_file(path_or_text)
    return parse_gene_codes_text(path_or_text)


def load_default_gene_codes() -> Dict[int, GeneCode]:
    """Load the bundled NCBI gc.prt resource."""
    try:
        ref = _ir.files('abaqs.data').joinpath('gc.prt')
        text = ref.read_text(encoding='utf-8')
    except Exception:
        _data_dir = Path(__file__).parent / 'data' / 'gc.prt'
        text = _data_dir.read_text(encoding='utf-8')
    return parse_gene_codes_text(text)


# =============================================================================
# Section 6 — Gff3Record
# =============================================================================

class Gff3Record:
    """Single GFF3 feature record with parent-child hierarchy support.

    Mirrors org.mycocosm.gff3.Gff3Record.
    """

    ATTRIBUTE_ID           = 'ID'
    ATTRIBUTE_PARENT       = 'Parent'
    ATTRIBUTE_PROTEIN_ID   = 'proteinId'
    ATTRIBUTE_TRANSCRIPT_ID = 'transcriptId'

    def __init__(
        self,
        category:  Gff3RecordCategory,
        id_:       Optional[str],
        seqid:     Optional[str],
        source:    Optional[str],
        type_:     Optional[Gff3Type],
        start:     int,
        end:       int,
        score:     float,
        strand:    Optional[Strand],
        phase:     int,
        attributes: Optional[Dict[str, str]] = None,
        meta_key:  Optional[str] = None,
        meta_val:  Optional[str] = None,
    ):
        self.category   = category
        self.id         = id_
        self.seqid      = seqid
        self.source     = source
        self.type       = type_
        self.start      = start
        self.end        = end
        self.score      = score
        self.strand     = strand
        self.phase      = phase
        self.attributes: Dict[str, str] = attributes or {}
        self.meta_key   = meta_key
        self.meta_val   = meta_val
        self.parent:   Optional['Gff3Record']  = None
        self.children: List['Gff3Record']      = []

    # ---- Factory helpers --------------------------------------------------

    @staticmethod
    def _regular(id_, seqid, source, type_, start, end, score, strand, phase, attributes):
        return Gff3Record(Gff3RecordCategory.regular, id_, seqid, source, type_,
                          start, end, score, strand, phase, attributes)

    @staticmethod
    def _meta(key: str, value: str) -> 'Gff3Record':
        return Gff3Record(Gff3RecordCategory.meta, None, None, None, None,
                          0, 0, float('nan'), None, -1,
                          meta_key=key, meta_val=value)

    @staticmethod
    def _fasta_marker() -> 'Gff3Record':
        return Gff3Record(Gff3RecordCategory.fasta, None, None, None, None,
                          0, 0, float('nan'), None, -1)

    # ---- Parsing ---------------------------------------------------------

    @staticmethod
    def from_line(line: str, accept_missing_id: bool = False) -> Optional['Gff3Record']:
        # Mirrors Java Gff3Record.fromLine — wraps every per-line failure with
        # the original line as context (A5), matching ExceptionsHelper.
        line = line.rstrip('\n\r')
        if not line or _COMMENT_RE.fullmatch(line):
            return None
        try:
            m = _META_RE.fullmatch(line)
            if m:
                key = _unescape(m.group(1))
                val = _unescape(m.group(2))
                if key == 'FASTA':
                    return Gff3Record._fasta_marker()
                return Gff3Record._meta(key, val)
            parts = line.split('\t')
            if len(parts) < 8:
                return None
            seqid   = _unescape(parts[0])
            source  = _unescape(parts[1])
            type_   = Gff3Type.value_of_or_none(parts[2])
            start   = _parse_position(parts[3])
            end     = _parse_position(parts[4])
            score   = _parse_score(parts[5])
            strand  = Strand.of(parts[6])
            phase   = _parse_phase(parts[7])
            attrs   = _parse_attributes(parts[8]) if len(parts) > 8 else {}
            id_     = attrs.get(Gff3Record.ATTRIBUTE_ID)
            if id_ is not None:
                return Gff3Record._regular(id_, seqid, source, type_,
                                           start, end, score, strand, phase, attrs)
            elif accept_missing_id:
                return Gff3Record._regular(None, seqid, source, type_,
                                           start, end, score, strand, phase, attrs)
            else:
                raise ValueError(f"Missing ID attribute in GFF3 line: {line!r}")
        except (ValueError, IndexError) as e:
            raise ValueError(
                f"Exception while parsing GFF3 line: {line!r} -- {e}"
            ) from e

    # ---- Parent-child management -----------------------------------------

    def _set_parent(self, parent: 'Gff3Record') -> None:
        self.parent = parent
        parent.children.append(self)

    # ---- Recursive search ------------------------------------------------

    def get_all_by_predicate(self,
                              predicate: Callable[['Gff3Record'], bool]
                              ) -> List['Gff3Record']:
        result: List['Gff3Record'] = []
        self._collect(predicate, result)
        return result

    def _collect(self, predicate, result):
        if predicate(self):
            result.append(self)
        for child in self.children:
            child._collect(predicate, result)

    @staticmethod
    def filter_top_level_regular(rec: 'Gff3Record') -> bool:
        return rec.category == Gff3RecordCategory.regular and rec.parent is None

    # ---- Genomic-sequence helpers ----------------------------------------

    def get_segments_by_type(self, type_: Gff3Type) -> List[Tuple[int, int]]:
        """Return ``(start, end)`` segments of the given type, sorted by start."""
        recs = self.get_all_by_predicate(lambda r: r.type == type_)
        recs.sort(key=lambda r: (r.start, r.end))
        return [(r.start, r.end) for r in recs]

    def get_translated_aminoacid_sequence(self, scaffold_seq: str, code: GeneCode
                                          ) -> Tuple[str, Optional[str]]:
        """Translate CDS to protein. Returns ``(protein_seq, overhang)``.

        Mirrors ``Gff3Record.getTranslatedAminoacidSequence()``. CDS segments
        are sorted by position, concatenated, then reverse-complemented if on
        the minus strand, then translated.
        """
        cds_segs = self.get_segments_by_type(Gff3Type.CDS)
        cds = ''.join(scaffold_seq[s - 1:e] for s, e in cds_segs)
        if self.strand == Strand.minus:
            cds = _reverse_complement(cds)
        return translate_sequence(cds, code)

    def get_transcript_sequence(self, scaffold_seq: str) -> str:
        """Return cDNA sequence (exon concatenation, oriented)."""
        exon_segs = self.get_segments_by_type(Gff3Type.exon)
        seq = ''.join(scaffold_seq[s - 1:e] for s, e in exon_segs)
        if self.strand == Strand.minus:
            seq = _reverse_complement(seq)
        return seq

    def length(self) -> int:
        return abs(self.end - self.start)

    # ---- GFF3 output formatting ------------------------------------------

    def _format_score(self) -> str:
        if math.isnan(self.score):
            return '.'
        if self.score == 0.0:
            return '0'
        return str(self.score)

    def _format_phase(self) -> str:
        return '.' if self.phase < 0 else str(self.phase)

    def _format_attributes(self) -> str:
        parts: List[str] = []
        if self.id is not None:
            parts.append(f"{_escape(self.ATTRIBUTE_ID)}={_escape(self.id)}")
        if self.parent is not None:
            parts.append(f"{_escape(self.ATTRIBUTE_PARENT)}={_escape(self.parent.id)}")
        for k, v in sorted(self.attributes.items(), key=lambda e: e[0].lower()):
            if k not in (self.ATTRIBUTE_ID, self.ATTRIBUTE_PARENT):
                parts.append(f"{_escape(k)}={_escape(v)}")
        return ';'.join(parts)

    def write(self, out) -> None:
        """Write this record (and children) in GFF3 format."""
        if self.category == Gff3RecordCategory.regular:
            out.write(
                f"{_escape(self.seqid)}\t{_escape(self.source)}\t"
                f"{self.type.value if self.type else '.'}\t"
                f"{self.start}\t{self.end}\t"
                f"{self._format_score()}\t"
                f"{self.strand.value if self.strand else '.'}\t"
                f"{self._format_phase()}\t"
                f"{self._format_attributes()}\n"
            )
            for child in self.children:
                child.write(out)
        elif self.category == Gff3RecordCategory.meta:
            if self.meta_val:
                out.write(f"##{_escape(self.meta_key)} {_escape(self.meta_val)}\n")
            else:
                out.write(f"##{_escape(self.meta_key)}\n")
        elif self.category == Gff3RecordCategory.fasta:
            out.write("##FASTA\n")


# =============================================================================
# Section 7 — GFF3Data
# =============================================================================

class GFF3Data:
    """Container for parsed GFF3 data.

    Mirrors org.mycocosm.gff3.GFF3Data. ``records`` contains only top-level
    records (no parent). ``scaffolds`` maps scaffold name → DNA sequence string.
    """

    def __init__(self,
                 records: List[Gff3Record],
                 scaffolds: Optional[Dict[str, str]]):
        self.records  = records
        self.scaffolds: Dict[str, str] = scaffolds or {}

    @staticmethod
    def parse(reader) -> 'GFF3Data':
        """Parse GFF3 from a text-mode iterable (file or StringIO)."""
        all_records: List[Gff3Record] = []
        by_id: Dict[str, Gff3Record]  = {}
        scaffolds: Optional[Dict[str, str]] = None

        for line in reader:
            if isinstance(line, bytes):
                line = line.decode()
            rec = Gff3Record.from_line(line)
            if rec is None:
                continue
            all_records.append(rec)
            if rec.id:
                by_id[rec.id] = rec
            if rec.category == Gff3RecordCategory.fasta:
                scaffolds = _load_fasta_from_reader(reader)
                break

        # Second pass: wire parent-child relationships
        top_level: List[Gff3Record] = []
        for rec in all_records:
            if rec.category != Gff3RecordCategory.regular:
                top_level.append(rec)
                continue
            parent_id = rec.attributes.get(Gff3Record.ATTRIBUTE_PARENT)
            if parent_id is None:
                top_level.append(rec)
            else:
                parent = by_id.get(parent_id)
                if parent is not None and parent is not rec:
                    rec._set_parent(parent)
                else:
                    top_level.append(rec)  # orphan → promote to top-level

        return GFF3Data(top_level, scaffolds)

    def has_scaffolds(self) -> bool:
        return bool(self.scaffolds)

    def get_records_by_predicate(self,
                                  predicate: Callable[[Gff3Record], bool]
                                  ) -> List[Gff3Record]:
        """Recursively collect all records (including children) matching predicate."""
        result: List[Gff3Record] = []
        for rec in self.records:
            result.extend(rec.get_all_by_predicate(predicate))
        return result

    def replace_scaffolds(self, new_scaffolds: Dict[str, str]) -> 'GFF3Data':
        return GFF3Data(self.records, new_scaffolds)

    def write_all(self, out, fasta_width: int = 70) -> None:
        for rec in self.records:
            rec.write(out)
            if rec.category == Gff3RecordCategory.fasta:
                for name, seq in self.scaffolds.items():
                    out.write(f">{name}\n")
                    for i in range(0, len(seq), fasta_width):
                        out.write(seq[i:i + fasta_width] + '\n')


# =============================================================================
# Section 8 — FASTA loading + IO helpers
# =============================================================================

def _load_fasta_from_reader(reader) -> Dict[str, str]:
    """Read FASTA from the current position in *reader*.

    Mirrors GFF3Data.loadFastaFromReader().
    """
    scaffolds: Dict[str, str] = {}
    current_name: Optional[str] = None
    current_seq: List[str] = []

    for line in reader:
        if isinstance(line, bytes):
            line = line.decode()
        line = line.rstrip('\n\r')
        if line.startswith('>') or line.startswith(';'):
            if current_name is not None:
                scaffolds[current_name] = ''.join(current_seq)
            header_parts = line[1:].split(None, 1)
            current_name = header_parts[0] if header_parts else ''
            current_seq = []
        else:
            current_seq.append(line)

    if current_name is not None:
        scaffolds[current_name] = ''.join(current_seq)

    return scaffolds


def load_fasta(path: str) -> Dict[str, str]:
    """Load a FASTA file (plain or gzipped) into a name→sequence dict."""
    opener = gzip.open if path.endswith('.gz') else open
    with opener(path, 'rt') as fh:
        return _load_fasta_from_reader(fh)


def open_gzipped_or_plain(path: str):
    """Return a text-mode reader for a possibly gzipped file."""
    if path.endswith('.gz'):
        return gzip.open(path, 'rt')
    return open(path, 'r')


def parse_gff3_file(path: str) -> 'GFF3Data':
    """Convenience: open (optionally gzipped) GFF3 file and parse it."""
    with open_gzipped_or_plain(path) as fh:
        return GFF3Data.parse(fh)


# =============================================================================
# Section 9 — Fixed-bin-size collection binning
# =============================================================================

T = TypeVar('T')


def bin_collection_fixed_bin_size(
    data: List[T],
    to_bin_value: Callable[[T], float],
    bin_step: float,
    min_value: Optional[float] = None,
    max_value: Optional[float] = None,
) -> List[Tuple[Tuple[float, float], List[T]]]:
    """Bin *data* into fixed-size bins.

    Returns a list of ``((from_bin_value, to_bin_value), [items])`` tuples,
    one per bin (including empty bins).

    ``min_value=None`` → use the minimum item value (Java Double.NaN → data min).
    ``max_value=None`` → use the maximum item value (Java Double.NaN → data max).

    Mirrors CollectionsHelper.binCollectionWithFixedBinSize().
    """
    if not data:
        return []

    entries = sorted(((to_bin_value(item), item) for item in data),
                     key=lambda x: x[0])

    min_bin = min_value if min_value is not None else entries[0][0]
    max_bin = max_value if max_value is not None else entries[-1][0]

    num_bins = int(math.floor((max_bin - min_bin) / bin_step))

    bins: List[List[T]] = [[] for _ in range(num_bins + 1)]
    ranges: List[Tuple[float, float]] = [
        (min_bin + i * bin_step, min_bin + (i + 1) * bin_step)
        for i in range(num_bins + 1)
    ]

    for val, item in entries:
        bin_idx = int(math.floor((val - min_bin) / bin_step))
        if 0 <= bin_idx < len(bins):
            bins[bin_idx].append(item)

    return list(zip(ranges, bins))
