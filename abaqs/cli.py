"""Command-line interface for ABAQS.

Mirrors the Java CLI flags defined in org.mycocosm.abaqs.main.ABAQS.
All short flags (-ig, -is, -ip, …) are preserved for drop-in compatibility.
"""

from __future__ import annotations

import argparse
import logging
import sys
from typing import Optional

from abaqs.scoring import (
    DEFAULT_BUSCO_LINEAGE,
    DEFAULT_BUSCO_THREADS,
    DEFAULT_GENE_CODE_ID,
    DEFAULT_ISOFORMS_MIN_OVERLAP,
    DEFAULT_MASKER_FUNCTION,
    DEFAULT_NO_DOMAINS_CDS_MASKED_CUTOFF,
    DEFAULT_PROTEIN_LENGTH_BINNING,
    DEFAULT_SUSPECTED_DOMAINS_MASKED_CUTOFF,
    compute_abaqs,
    print_result,
)

_VERSION = '1.0'


def _build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog='abaqs',
        description='Accuracy-Based Annotation Quality Score (ABAQS) — Python port',
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument('--version', action='version', version=f'ABAQS {_VERSION}')

    # ---- Required ----
    p.add_argument('-ig', '--input-gff', required=True, metavar='FILE',
                   help='Input GFF3 file path (plain or .gz)')

    # ---- Optional inputs ----
    p.add_argument('-is', '--input-scaffolds-fasta', metavar='FILE',
                   help='Input scaffolds FASTA (softmasked, plain or .gz)')
    p.add_argument('-ip', '--input-proteins-fasta', metavar='FILE',
                   help='Input proteins FASTA (plain or .gz); if omitted, CDS are translated')
    p.add_argument('-id', '--input-domains', metavar='FILE',
                   help='Input InterPro/Pfam domains file (plain or .gz)')
    p.add_argument('-ibf', '--busco-data-file', metavar='FILE',
                   help=("BUSCO short_summary.txt OR compleasm summary.txt file. "
                         "Format is auto-detected from the filename."))
    p.add_argument('-ib', '--busco-data', metavar='STRING',
                   help="BUSCO summary string, e.g. 'C:99.3%%[S:98.9%%,D:0.4%%],F:0.3%%,M:0.4%%,n:758'")
    p.add_argument('--no-busco-auto-run', dest='busco_auto_run',
                   action='store_false', default=True,
                   help=("Disable auto-discovery and auto-run of compleasm/BUSCO when "
                         "neither --busco-data nor --busco-data-file is provided. "
                         "By default, ABAQS looks for "
                         "'<basename>_<lineage>_compleasm/summary.txt' beside "
                         "--input-scaffolds-fasta, falling back to "
                         "'<basename>_<lineage>_busco/short_summary*.txt', and then "
                         "to invoking compleasm/busco as a subprocess."))
    p.add_argument('--busco-lineage', metavar='LINEAGE',
                   default=DEFAULT_BUSCO_LINEAGE,
                   help=(f"Lineage to prefer when discovering or running compleasm/BUSCO "
                         f"(default: {DEFAULT_BUSCO_LINEAGE})."))
    p.add_argument('--busco-threads', type=int, metavar='INT',
                   default=DEFAULT_BUSCO_THREADS,
                   help=f"Threads to pass to compleasm/BUSCO (default: {DEFAULT_BUSCO_THREADS}).")
    p.add_argument('-ite', '--te-input-file', metavar='FILE',
                   help='Transposable-element Pfam domain list (uses bundled list if omitted)')
    p.add_argument('-ise', '--suspected-te-input-file', metavar='FILE',
                   help='Suspected TE Pfam domain list (uses bundled list if omitted)')
    p.add_argument('-ilr', '--reference-protein-lengths-input-file', metavar='FILE',
                   help='Reference protein-length distribution (uses bundled data if omitted)')
    p.add_argument('-igc', '--gene-code-input-file', metavar='FILE',
                   help='NCBI gc.prt gene-code file (uses bundled file if omitted)')

    # ---- Optional outputs ----
    p.add_argument('-o', '--output', metavar='FILE',
                   help='Output path for results (default: stdout)')
    p.add_argument('-og', '--output-gff', metavar='FILE',
                   help='Write parsed GFF3 to this file')

    # ---- Mapper patterns ----
    p.add_argument('-mg', '--gff3-protein-id-mapper',
                   default='proteinId', metavar='ATTR_OR_SPEC',
                   help=("GFF3 protein-ID mapper. Either an attribute name (e.g. 'proteinId') "
                         "or the Java form 'column:attribute:pattern->template' where column "
                         "is seqid/source/attributes (default: 'proteinId')"))
    p.add_argument('-mp', '--protein-fasta-protein-id-mapper', metavar='PATTERN',
                   help=(r"PatternTransformer for protein FASTA ID, format 'regex->template' "
                         r"with optional '<<CI,DA' flags; supports {0..n} and {name} groups "
                         r"(default: .+proteinId\s*=\s*(\d+).*->{1})"))
    p.add_argument('-md', '--domains-protein-id-mapper', metavar='PATTERN',
                   help='Regex with named groups <id> and <domain> for the domain TSV (default: HMMPfam TSV)')

    # ---- Numeric tuning ----
    p.add_argument('-g', '--gene-code', type=int, default=DEFAULT_GENE_CODE_ID,
                   metavar='INT',
                   help=f'NCBI genetic code ID (default: {DEFAULT_GENE_CODE_ID})')
    p.add_argument('-io', '--isoforms-min-overlap', type=float,
                   default=DEFAULT_ISOFORMS_MIN_OVERLAP, metavar='FLOAT',
                   help=f'Min overlap to detect isoforms (default: {DEFAULT_ISOFORMS_MIN_OVERLAP})')
    p.add_argument('-mf', '--masker-function',
                   default=DEFAULT_MASKER_FUNCTION,
                   choices=['TO_LOWER_CASE', 'TO_N', 'TO_DOT'],
                   help=f'Soft-masking detection function (default: {DEFAULT_MASKER_FUNCTION})')
    p.add_argument('-ndc', '--no-domain-masked-cutoff', type=float,
                   default=DEFAULT_NO_DOMAINS_CDS_MASKED_CUTOFF, metavar='FLOAT',
                   help=f'CDS-masked cutoff for no-domain TE detection (default: {DEFAULT_NO_DOMAINS_CDS_MASKED_CUTOFF}; NaN = disabled)')
    p.add_argument('-sdc', '--suspected-domain-masked-cutoff', type=float,
                   default=DEFAULT_SUSPECTED_DOMAINS_MASKED_CUTOFF, metavar='FLOAT',
                   help=f'CDS-masked cutoff for suspected-TE detection (default: NaN = always TE)')
    p.add_argument('-plb', '--protein-length-binning', type=int,
                   default=DEFAULT_PROTEIN_LENGTH_BINNING, metavar='INT',
                   help=f'Protein length bin size (default: {DEFAULT_PROTEIN_LENGTH_BINNING})')
    p.add_argument('-fw', '--fasta-width', type=int, default=70, metavar='INT',
                   help='FASTA output line width (default: 70)')

    # ---- Verbosity ----
    p.add_argument('-v', '--verbose', action='store_true',
                   help='Verbose logging output')
    p.add_argument('-vo', '--verbose-output-folder', metavar='DIR',
                   help='Write verbose computation files to this folder')

    return p


def main(argv=None) -> None:
    parser = _build_parser()
    args = parser.parse_args(argv)

    log_level = logging.DEBUG if args.verbose else logging.INFO
    logging.basicConfig(
        level=log_level,
        format='%(asctime)s %(levelname)s %(name)s: %(message)s',
        stream=sys.stderr,
    )

    # ---- Run ----
    result = compute_abaqs(
        input_gff3                    = args.input_gff,
        input_scaffolds_fasta         = args.input_scaffolds_fasta,
        input_proteins_fasta          = args.input_proteins_fasta,
        input_domains                 = args.input_domains,
        busco_string                  = args.busco_data,
        busco_file                    = args.busco_data_file,
        te_domains_file               = args.te_input_file,
        suspected_te_domains_file     = args.suspected_te_input_file,
        reference_pld_file            = args.reference_protein_lengths_input_file,
        gene_code_file                = args.gene_code_input_file,
        gene_code_id                  = args.gene_code,
        gff3_protein_id_attr          = args.gff3_protein_id_mapper,
        protein_fasta_id_pattern      = args.protein_fasta_protein_id_mapper,
        domains_pattern               = args.domains_protein_id_mapper,
        isoforms_min_overlap          = args.isoforms_min_overlap,
        masker_function               = args.masker_function,
        no_domain_cds_masked_cutoff   = args.no_domain_masked_cutoff,
        suspected_domain_masked_cutoff = args.suspected_domain_masked_cutoff,
        protein_length_binning        = args.protein_length_binning,
        busco_auto_run                = args.busco_auto_run,
        busco_lineage                 = args.busco_lineage,
        busco_threads                 = args.busco_threads,
        verbose                       = args.verbose,
    )

    # ---- Output GFF3 ----
    if args.output_gff:
        import gzip as _gz
        from abaqs.sequence import parse_gff3_file
        gff3_out = parse_gff3_file(args.input_gff)
        opener = _gz.open if args.output_gff.endswith('.gz') else open
        with opener(args.output_gff, 'wt') as fh:
            gff3_out.write_all(fh, fasta_width=args.fasta_width)

    # ---- Print result ----
    if args.output:
        import gzip as _gz
        opener = _gz.open if args.output.endswith('.gz') else open
        with opener(args.output, 'wt') as fh:
            print_result(result, args.input_gff, fh)
    else:
        print_result(result, args.input_gff)
