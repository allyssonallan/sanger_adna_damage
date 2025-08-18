"""
Primer detection and handling components for Enhanced AB1 Converter.

This module provides specialized functionality for detecting primer orientations,
removing primers, and handling F/R primer pairs.
"""

import logging
from typing import Dict, Any, Tuple, Optional

logger = logging.getLogger(__name__)


class PrimerHandler:
    """Handles primer detection, orientation, and removal operations."""

    def __init__(
        self,
        custom_primers_forward: Optional[Dict[str, str]] = None,
        custom_primers_reverse: Optional[Dict[str, str]] = None,
    ):
        """
        Initialize primer handler.

        Args:
            custom_primers_forward: Optional custom forward primers
            custom_primers_reverse: Optional custom reverse primers
        """
        self.primer_pairs = self._setup_primer_pairs(
            custom_primers_forward, custom_primers_reverse
        )

    def _setup_primer_pairs(
        self,
        custom_forward: Optional[Dict[str, str]] = None,
        custom_reverse: Optional[Dict[str, str]] = None,
    ) -> Dict[str, Dict[str, str]]:
        """
        Set up primer pairs for HVS regions with F/R capability.

        Args:
            custom_forward: Optional custom forward primers
            custom_reverse: Optional custom reverse primers

        Returns:
            Dictionary of primer pairs organized by HVS region
        """
        # Default primer sequences (these are examples - should be updated with actual sequences)
        default_primers = {
            "HVS1": {
                "forward": "TTCTTTCATGGGGAAGCAGATTTGGGT",
                "reverse": "TGATTTTTTGGATGTAGATGAGGAG",
            },
            "HVS2": {
                "forward": "CTGACCGTGAAATCAGCAAAC",
                "reverse": "AGTGGATGAGGCAGGTTGTT",
            },
            "HVS3": {
                "forward": "CACCATGAATATCATTGGTC",
                "reverse": "GGTGATGTGAAGGTTGATGT",
            },
        }

        # Merge with custom primers if provided
        primers = default_primers.copy()

        if custom_forward:
            for region, forward_primer in custom_forward.items():
                if region in primers:
                    primers[region]["forward"] = forward_primer
                else:
                    primers[region] = {"forward": forward_primer, "reverse": ""}

        if custom_reverse:
            for region, reverse_primer in custom_reverse.items():
                if region in primers:
                    primers[region]["reverse"] = reverse_primer
                else:
                    primers[region] = {"forward": "", "reverse": reverse_primer}

        return primers

    def detect_primer_orientation(self, sequence: str, region: str) -> Dict[str, Any]:
        """
        Detect primer orientation (F or R) in sequence.

        Args:
            sequence: DNA sequence to analyze
            region: HVS region (HVS1, HVS2, HVS3)

        Returns:
            Dictionary with orientation detection results
        """
        if region not in self.primer_pairs:
            return {
                "orientation": "unknown",
                "confidence": 0.0,
                "primer_found": False,
                "details": f"No primers defined for region {region}",
            }

        primers = self.primer_pairs[region]
        forward_primer = primers.get("forward", "")
        reverse_primer = primers.get("reverse", "")

        # Normalize sequence
        norm_sequence = self._normalize_sequence(sequence)

        # Check for forward primer
        forward_score = 0.0
        if forward_primer:
            forward_score = self._find_best_primer_match(
                norm_sequence, forward_primer, from_start=True
            )

        # Check for reverse primer (and its reverse complement)
        reverse_score = 0.0
        if reverse_primer:
            reverse_score = self._find_best_primer_match(
                norm_sequence, reverse_primer, from_start=True
            )
            # Also check reverse complement of reverse primer
            reverse_rc = self._reverse_complement(reverse_primer)
            reverse_rc_score = self._find_best_primer_match(
                norm_sequence, reverse_rc, from_start=True
            )
            reverse_score = max(reverse_score, reverse_rc_score)

        # Determine orientation
        if forward_score > reverse_score and forward_score > 0.7:
            orientation = "F"
            confidence = forward_score
        elif reverse_score > forward_score and reverse_score > 0.7:
            orientation = "R"
            confidence = reverse_score
        else:
            orientation = "unknown"
            confidence = max(forward_score, reverse_score)

        return {
            "orientation": orientation,
            "confidence": confidence,
            "primer_found": confidence > 0.7,
            "details": {
                "forward_score": forward_score,
                "reverse_score": reverse_score,
                "region": region,
            },
        }

    def remove_primers_with_orientation(
        self, sequence: str, region: str
    ) -> Tuple[str, Dict]:
        """
        Remove primers based on detected orientation.

        Args:
            sequence: DNA sequence
            region: HVS region

        Returns:
            Tuple of (cleaned_sequence, removal_stats)
        """
        # Detect orientation first
        orientation_result = self.detect_primer_orientation(sequence, region)

        if not orientation_result["primer_found"]:
            return sequence, {
                "primers_removed": False,
                "orientation": "unknown",
                "original_length": len(sequence),
                "final_length": len(sequence),
                "removed_bases": 0,
            }

        orientation = orientation_result["orientation"]
        primers = self.primer_pairs.get(region, {})

        cleaned_sequence = sequence
        total_removed = 0

        if orientation == "F" and primers.get("forward"):
            # Remove forward primer from start
            primer_pos = self._find_primer_position(
                cleaned_sequence, primers["forward"], from_start=True
            )
            if primer_pos >= 0:
                primer_end = primer_pos + len(primers["forward"])
                cleaned_sequence = cleaned_sequence[primer_end:]
                total_removed += primer_end

        elif orientation == "R" and primers.get("reverse"):
            # For reverse orientation, we need to handle reverse complement
            # First, reverse complement the sequence
            cleaned_sequence = self._reverse_complement(cleaned_sequence)

            # Then remove the reverse primer (which should now be at the start)
            primer_pos = self._find_primer_position(
                cleaned_sequence, primers["reverse"], from_start=True
            )
            if primer_pos >= 0:
                primer_end = primer_pos + len(primers["reverse"])
                cleaned_sequence = cleaned_sequence[primer_end:]
                total_removed += primer_end

        return cleaned_sequence, {
            "primers_removed": total_removed > 0,
            "orientation": orientation,
            "original_length": len(sequence),
            "final_length": len(cleaned_sequence),
            "removed_bases": total_removed,
            "confidence": orientation_result["confidence"],
        }

    def _normalize_sequence(self, sequence: str) -> str:
        """
        Normalize sequence for primer matching.

        Args:
            sequence: Raw DNA sequence

        Returns:
            Normalized sequence
        """
        # Convert to uppercase and remove non-ATCG characters except N
        normalized = "".join(
            char.upper() for char in sequence if char.upper() in "ATCGN"
        )

        # Replace runs of N with single N for better matching
        import re

        normalized = re.sub(r"N+", "N", normalized)

        return normalized

    def _find_primer_position(
        self, sequence: str, primer: str, from_start: bool = True
    ) -> int:
        """
        Find primer position in sequence with fuzzy matching.

        Args:
            sequence: Target sequence
            primer: Primer sequence to find
            from_start: Whether to search from start (True) or end (False)

        Returns:
            Position of primer (-1 if not found)
        """
        if not primer:
            return -1

        norm_sequence = self._normalize_sequence(sequence)
        norm_primer = self._normalize_sequence(primer)

        # Simple exact match first
        if from_start:
            pos = norm_sequence.find(norm_primer)
        else:
            pos = norm_sequence.rfind(norm_primer)

        if pos >= 0:
            return pos

        # Fuzzy matching with up to 2 mismatches
        min_length = min(len(norm_sequence), len(norm_primer))
        search_length = min(50, min_length)  # Limit search window

        best_score = 0.0
        best_pos = -1

        if from_start:
            search_range = range(
                0, min(search_length, len(norm_sequence) - len(norm_primer) + 1)
            )
        else:
            start_pos = max(0, len(norm_sequence) - search_length)
            search_range = range(start_pos, len(norm_sequence) - len(norm_primer) + 1)

        for pos in search_range:
            if pos + len(norm_primer) <= len(norm_sequence):
                subseq = norm_sequence[pos : pos + len(norm_primer)]
                score = self._calculate_match_score(subseq, norm_primer)

                if score > best_score and score > 0.8:  # 80% similarity threshold
                    best_score = score
                    best_pos = pos

        return best_pos if best_score > 0.8 else -1

    def _find_best_primer_match(
        self, sequence: str, primer: str, from_start: bool = True
    ) -> float:
        """
        Find best primer match score in sequence.

        Args:
            sequence: Target sequence
            primer: Primer sequence to match
            from_start: Whether to search from start

        Returns:
            Best match score (0.0 to 1.0)
        """
        pos = self._find_primer_position(sequence, primer, from_start)
        if pos < 0:
            return 0.0

        norm_sequence = self._normalize_sequence(sequence)
        norm_primer = self._normalize_sequence(primer)

        if pos + len(norm_primer) <= len(norm_sequence):
            subseq = norm_sequence[pos : pos + len(norm_primer)]
            return self._calculate_match_score(subseq, norm_primer)

        return 0.0

    def _calculate_match_score(self, seq1: str, seq2: str) -> float:
        """
        Calculate similarity score between two sequences.

        Args:
            seq1: First sequence
            seq2: Second sequence

        Returns:
            Similarity score (0.0 to 1.0)
        """
        if not seq1 or not seq2:
            return 0.0

        min_len = min(len(seq1), len(seq2))
        if min_len == 0:
            return 0.0

        matches = sum(1 for i in range(min_len) if seq1[i] == seq2[i])
        return matches / min_len

    def _reverse_complement(self, sequence: str) -> str:
        """
        Generate reverse complement of DNA sequence.

        Args:
            sequence: DNA sequence

        Returns:
            Reverse complement sequence
        """
        complement_map = {
            "A": "T",
            "T": "A",
            "C": "G",
            "G": "C",
            "N": "N",
            "R": "Y",
            "Y": "R",
            "M": "K",
            "K": "M",
            "S": "S",
            "W": "W",
            "B": "V",
            "V": "B",
            "D": "H",
            "H": "D",
        }

        # Convert to uppercase and get complement
        complement = "".join(
            complement_map.get(base.upper(), base) for base in sequence
        )

        # Return reverse complement
        return complement[::-1]
