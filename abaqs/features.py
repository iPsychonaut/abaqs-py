"""Genomic feature representation and track-to-track overlap mapping.

Mirrors:
  org.mycocosm.framework.biology.{FeatureGeneric,FeatureModel,Locus,LocusKey}
  org.mycocosm.abaqs.compare_tracks.{FeaturePair,FeatureTrack,TrackToTrackMapping}

Intentional divergences from Java (Python is more correct):
  - B2: A model feature with no CDS falls back to the mRNA genomic span
        instead of Integer.MAX_VALUE / Integer.MIN_VALUE sentinels.
  - B3: ``TrackToTrackMapping.percent_mapped`` uses the same denominator as
        Java (``mapped + side_a_not + side_b_not``) — earlier port omitted
        ``side_b_not_mapped``.
  - C1: ``map_two_tracks_by_position_overlap`` iterates ``track_a.features``
        in deterministic order (sorted by ``feature_id``). Java's HashSet
        iteration is non-deterministic; tie resolution during displacement
        contests previously depended on hash order.
"""

from __future__ import annotations

import functools
import logging
from dataclasses import dataclass, field
from typing import Callable, Dict, List, Optional, Set, Tuple

from abaqs.sequence import GFF3Data, Gff3Record, Gff3Type, Strand

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# LocusKey  (scaffold + strand, used as dict key)
# ---------------------------------------------------------------------------

class LocusKey:
    """Immutable key combining scaffold name and strand.

    Mirrors org.mycocosm.framework.biology.LocusKey.
    """

    __slots__ = ('scaffold', 'strand', '_hash')

    def __init__(self, scaffold: str, strand: Strand):
        self.scaffold = scaffold
        self.strand   = strand
        self._hash    = hash((scaffold, strand))

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, LocusKey):
            return NotImplemented
        return self.scaffold == other.scaffold and self.strand == other.strand

    def __hash__(self) -> int:
        return self._hash

    def __repr__(self) -> str:
        return f"LocusKey({self.scaffold!r}, {self.strand})"


# ---------------------------------------------------------------------------
# Feature  (combines FeatureGeneric + FeatureModel in a single class)
# ---------------------------------------------------------------------------

class Feature:
    """A genomic feature built from a GFF3 mRNA record.

    In *generic* mode  effective_start/end = genomic span (mRNA coords).
    In *model*   mode  effective_start/end = CDS span.

    Mirrors FeatureGeneric / FeatureModel.
    """

    def __init__(
        self,
        feature_id:      int,          # internal numeric id (from counter)
        name:            str,          # mRNA.id
        locus_key:       LocusKey,
        genomic_start:   int,
        genomic_end:     int,
        effective_start: int,
        effective_end:   int,
        exon_starts:     List[int],
        exon_ends:       List[int],
    ):
        self.feature_id      = feature_id
        self.name            = name
        self.locus_key       = locus_key
        self.genomic_start   = genomic_start
        self.genomic_end     = genomic_end
        self.effective_start = effective_start
        self.effective_end   = effective_end
        self.exon_starts     = exon_starts
        self.exon_ends       = exon_ends

    # ------------------------------------------------------------------
    # Factory
    # ------------------------------------------------------------------

    @staticmethod
    def from_gff3_mrna_generic(mrna: Gff3Record, feature_id: int) -> 'Feature':
        """Build a feature using exon span as the effective coordinates."""
        exon_segs = mrna.get_segments_by_type(Gff3Type.exon)
        starts = [s for s, _ in exon_segs]
        ends   = [e for _, e in exon_segs]
        lkey   = LocusKey(mrna.seqid, mrna.strand)
        return Feature(
            feature_id      = feature_id,
            name            = mrna.id,
            locus_key       = lkey,
            genomic_start   = mrna.start,
            genomic_end     = mrna.end,
            effective_start = mrna.start,
            effective_end   = mrna.end,
            exon_starts     = starts,
            exon_ends       = ends,
        )

    @staticmethod
    def from_gff3_mrna_model(mrna: Gff3Record, feature_id: int) -> 'Feature':
        """Build a feature using CDS span as the effective coordinates."""
        exon_segs = mrna.get_segments_by_type(Gff3Type.exon)
        starts = [s for s, _ in exon_segs]
        ends   = [e for _, e in exon_segs]
        cds_segs = mrna.get_segments_by_type(Gff3Type.CDS)
        cds_start = min((s for s, _ in cds_segs), default=mrna.start)
        cds_end   = max((e for _, e in cds_segs), default=mrna.end)
        lkey = LocusKey(mrna.seqid, mrna.strand)
        return Feature(
            feature_id      = feature_id,
            name            = mrna.id,
            locus_key       = lkey,
            genomic_start   = mrna.start,
            genomic_end     = mrna.end,
            effective_start = cds_start,
            effective_end   = cds_end,
            exon_starts     = starts,
            exon_ends       = ends,
        )

    # ------------------------------------------------------------------
    # Overlap computation
    # ------------------------------------------------------------------

    def compute_overlap(self, other: 'Feature') -> float:
        """Normalized overlap in [0, 1].

        0.0  = no overlap; 1.0 = identical span.
        Only meaningful for same scaffold + strand.
        Mirrors FeatureGeneric.computeOverlap().
        """
        if (self.locus_key.scaffold != other.locus_key.scaffold or
                self.locus_key.strand != other.locus_key.strand):
            return 0.0
        a_from = min(self.effective_start, self.effective_end)
        a_to   = max(self.effective_start, self.effective_end)
        b_from = min(other.effective_start, other.effective_end)
        b_to   = max(other.effective_start, other.effective_end)
        if a_from < b_to and a_to > b_from:
            base    = max(a_to, b_to) - min(a_from, b_from)
            overlap = min(a_to, b_to) - max(a_from, b_from)
            return overlap / base
        return 0.0

    def __repr__(self) -> str:
        return (f"Feature(id={self.feature_id}, name={self.name!r}, "
                f"eff=[{self.effective_start},{self.effective_end}])")


# ---------------------------------------------------------------------------
# FeaturePair
# ---------------------------------------------------------------------------

@dataclass
class FeaturePair:
    """A matched pair (A, B) with their overlap ratio.

    Mirrors org.mycocosm.abaqs.compare_tracks.FeaturePair.
    """
    a: Feature
    b: Feature
    overlap: float

    @staticmethod
    def compare_by_overlap_desc(p1: 'FeaturePair', p2: 'FeaturePair') -> int:
        """Comparator: descending overlap (higher overlap = earlier in sort).

        Double.compare(f2.overlap, f1.overlap)
        """
        if p2.overlap < p1.overlap:
            return -1
        if p2.overlap > p1.overlap:
            return 1
        return 0


# ---------------------------------------------------------------------------
# TrackToTrackMapping
# ---------------------------------------------------------------------------

@dataclass
class TrackToTrackMapping:
    """Result of mapping two feature tracks by positional overlap.

    Mirrors org.mycocosm.abaqs.compare_tracks.TrackToTrackMapping.
    """
    mapped_pairs:    List[FeaturePair]
    side_a_not_mapped: Set[Feature]
    side_b_not_mapped: Set[Feature]

    def percent_mapped(self) -> float:
        # B3: Java denominator includes both unmapped sides.
        total = (len(self.mapped_pairs)
                 + len(self.side_a_not_mapped)
                 + len(self.side_b_not_mapped))
        return len(self.mapped_pairs) / total if total else 0.0


# ---------------------------------------------------------------------------
# FeatureTrack
# ---------------------------------------------------------------------------

class FeatureTrack:
    """A set of genomic features with a spatial index for overlap queries.

    Mirrors org.mycocosm.abaqs.compare_tracks.FeatureTrack.
    """

    def __init__(self, name: str, model_mode: bool = True):
        self.name       = name
        self.model_mode = model_mode          # True → use CDS coords
        self.features:  Set[Feature] = set()
        self._index:    Dict[LocusKey, List[Feature]] = {}  # sorted by eff_start

    # ------------------------------------------------------------------
    # Factory
    # ------------------------------------------------------------------

    @staticmethod
    def create_from_gff3(gff3: GFF3Data, name: str = 'track',
                         model: bool = True) -> 'FeatureTrack':
        """Build a track from all mRNA records in *gff3*.

        Mirrors FeatureTrack.createTrackFromGff3Data().
        """
        track    = FeatureTrack(name, model)
        counter  = 0
        for mrna in gff3.get_records_by_predicate(lambda r: r.type == Gff3Type.mRNA):
            counter += 1
            if model:
                feat = Feature.from_gff3_mrna_model(mrna, counter)
            else:
                feat = Feature.from_gff3_mrna_generic(mrna, counter)
            track.features.add(feat)
        track._build_index()
        return track

    def _build_index(self) -> None:
        """Create spatial index: LocusKey → features sorted by effective_start.

        Mirrors FeatureTrack.createFeaturesIndex().
        """
        from collections import defaultdict
        tmp: Dict[LocusKey, List[Feature]] = {}
        for feat in self.features:
            lk = feat.locus_key
            if lk not in tmp:
                tmp[lk] = []
            tmp[lk].append(feat)
        self._index = {lk: sorted(feats, key=lambda f: f.effective_start)
                       for lk, feats in tmp.items()}

    # ------------------------------------------------------------------
    # Track-to-track overlap mapping
    # ------------------------------------------------------------------

    @staticmethod
    def map_two_tracks_by_position_overlap(
        track_a: 'FeatureTrack',
        track_b: 'FeatureTrack',
        best_pair_comparator: Callable[[FeaturePair, FeaturePair], int],
    ) -> TrackToTrackMapping:
        """Bidirectional greedy overlap mapping with displacement.

        Mirrors FeatureTrack.mapTwoTracksByPositionOverlap().

        ``best_pair_comparator(p1, p2)`` should return < 0 when p1 is
        "preferred" over p2 (i.e. should survive a displacement contest).
        For isoform detection the caller passes compareByOverlapDesc so that
        lower-overlap pairs are "preferred" (more conservative detection).
        """
        mapped_pairs:          Dict[int, FeaturePair] = {}  # id(feat_b) → pair
        reconsider:            List[Feature]           = []
        # C1: deterministic iteration order (feature_id assigned by counter
        # while parsing the GFF3, so this reflects GFF3 source order).
        side_a_list:           List[Feature]           = sorted(
            track_a.features, key=lambda f: f.feature_id)

        while side_a_list:
            next_reconsider: List[Feature] = []
            for feat_a in side_a_list:
                pair = track_b._find_best_by_overlap(
                    feat_a, mapped_pairs, next_reconsider, best_pair_comparator
                )
                if pair is not None:
                    mapped_pairs[id(pair.b)] = pair
            side_a_list = next_reconsider

        all_pairs  = list(mapped_pairs.values())
        mapped_a   = {id(p.a) for p in all_pairs}
        mapped_b   = {id(p.b) for p in all_pairs}
        not_a      = {f for f in track_a.features if id(f) not in mapped_a}
        not_b      = {f for f in track_b.features if id(f) not in mapped_b}

        return TrackToTrackMapping(all_pairs, not_a, not_b)

    def _find_best_by_overlap(
        self,
        feat_a: Feature,
        mapped_pairs: Dict[int, FeaturePair],
        reconsider_queue: List[Feature],
        comparator: Callable[[FeaturePair, FeaturePair], int],
    ) -> Optional[FeaturePair]:
        """Find best matching feature in this track for *feat_a*.

        Mirrors FeatureTrack.findBestFeatureByOverlap().
        Uses the same binary-search + walk-back algorithm as the Java code.
        """
        segment = self._index.get(feat_a.locus_key)
        if not segment:
            return None
        n = len(segment)

        # ---- binary search for insertion point of feat_a.effective_end ----
        feature_end = feat_a.effective_end
        left  = 0
        right = n - 1
        mid   = left + (right - left) // 2

        found = False
        while not found:
            if left >= right - 1:
                found = True
                mid   = right
            else:
                mid_start  = segment[mid].effective_start
                down_start = segment[mid - 1].effective_start
                if down_start < feature_end and mid_start >= feature_end:
                    found = True
                elif down_start < feature_end and mid_start < feature_end:
                    left = mid
                    mid  = left + (right - left) // 2
                else:
                    right = mid
                    mid   = left + (right - left) // 2

        # ---- walk backward collecting all overlapping candidates ----
        current  = mid
        possible: List[FeaturePair] = []
        while current >= 0 and segment[current].effective_end > feat_a.effective_start:
            feat_b = segment[current]
            if feat_b is not feat_a:
                overlap = feat_b.compute_overlap(feat_a)
                if overlap > 0.0:
                    possible.append(FeaturePair(feat_a, feat_b, overlap))
            current -= 1

        if not possible:
            return None

        # Sort with REVERSED comparator (mirrors Java's ".sorted((p1,p2)->cmp.compare(p2,p1))")
        possible.sort(key=functools.cmp_to_key(lambda p1, p2: comparator(p2, p1)))

        for m in possible:
            prev = mapped_pairs.get(id(m.b))
            if prev is not None:
                if comparator(prev, m) < 0:
                    # m is "better" than prev → displace prev.a back to reconsider queue
                    del mapped_pairs[id(m.b)]
                    reconsider_queue.append(prev.a)
                    return m
                # else prev is better; try next candidate
            else:
                return m  # B is free → take it

        return None  # all candidates already mapped by better pairs
