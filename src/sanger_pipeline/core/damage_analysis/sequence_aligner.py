"""
Sequence alignment component for ancient DNA damage analysis.
"""

import logging
from typing import Tuple
from Bio.Align import PairwiseAligner

logger = logging.getLogger(__name__)


class SequenceAligner:
    """Handles sequence alignment for damage analysis."""

    def __init__(
        self,
        match_score: int = 2,
        mismatch_score: int = -1,
        open_gap_score: int = -2,
        extend_gap_score: float = -0.5,
    ):
        """
        Initialize sequence aligner.

        Args:
            match_score: Score for matching bases
            mismatch_score: Score for mismatching bases
            open_gap_score: Score for opening a gap
            extend_gap_score: Score for extending a gap
        """
        self.match_score = match_score
        self.mismatch_score = mismatch_score
        self.open_gap_score = open_gap_score
        self.extend_gap_score = extend_gap_score

    def align_sequences(self, reference: str, query: str) -> Tuple[str, str]:
        """
        Align two sequences using global alignment.

        Args:
            reference: Reference sequence string
            query: Query sequence string

        Returns:
            Tuple of (aligned_reference, aligned_query)
        """
        aligner = PairwiseAligner()
        aligner.match_score = self.match_score
        aligner.mismatch_score = self.mismatch_score
        aligner.open_gap_score = self.open_gap_score
        aligner.extend_gap_score = self.extend_gap_score
        aligner.mode = "global"

        alignments = aligner.align(reference, query)

        if not alignments:
            raise ValueError("No alignment found between sequences")

        # Take the best alignment (first one)
        best_alignment = alignments[0]
        aligned_ref = str(best_alignment[0])
        aligned_query = str(best_alignment[1])

        return aligned_ref, aligned_query

    def clean_alignment(self, aligned_ref: str, aligned_query: str) -> Tuple[str, str]:
        """
        Remove gaps from alignment, keeping only positions with reference bases.

        Args:
            aligned_ref: Aligned reference sequence
            aligned_query: Aligned query sequence

        Returns:
            Tuple of (clean_reference, clean_query)
        """
        ref_clean = ""
        query_clean = ""

        for ref_base, query_base in zip(aligned_ref, aligned_query):
            if ref_base != "-":  # Keep positions with reference bases
                ref_clean += ref_base
                query_clean += query_base if query_base != "-" else "N"

        return ref_clean, query_clean
