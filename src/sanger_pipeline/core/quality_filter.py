"""
Quality filtering utilities for sequence data.
"""

import logging
from typing import List

from ..utils.constants import DEFAULT_MIN_QUALITY


logger = logging.getLogger(__name__)


class QualityFilter:
    """
    Quality filtering for sequence data.
    """

    def __init__(self, min_quality: int = DEFAULT_MIN_QUALITY):
        """
        Initialize quality filter.

        Args:
            min_quality: Minimum quality score threshold
        """
        self.min_quality = min_quality

    def filter_sequence(self, sequence: str, qualities: List[int]) -> str:
        """
        Filter sequence based on quality scores.

        Args:
            sequence: DNA sequence string
            qualities: List of quality scores

        Returns:
            Filtered sequence with low-quality bases replaced by 'N'
        """
        if len(sequence) != len(qualities):
            raise ValueError("Sequence and qualities must have the same length")

        filtered_seq = "".join(
            [base if qual >= self.min_quality else "N" for base, qual in zip(sequence, qualities)]
        )

        return filtered_seq

    def calculate_quality_stats(self, qualities: List[int]) -> dict:
        """
        Calculate quality statistics.

        Args:
            qualities: List of quality scores

        Returns:
            Dictionary with quality statistics
        """
        if not qualities:
            return {}

        return {
            "mean_quality": sum(qualities) / len(qualities),
            "min_quality": min(qualities),
            "max_quality": max(qualities),
            "bases_above_threshold": sum(1 for q in qualities if q >= self.min_quality),
            "total_bases": len(qualities),
            "fraction_above_threshold": sum(1 for q in qualities if q >= self.min_quality)
            / len(qualities),
        }

    def get_quality_regions(self, qualities: List[int]) -> List[dict]:
        """
        Identify regions of high and low quality.

        Args:
            qualities: List of quality scores

        Returns:
            List of dictionaries describing quality regions
        """
        regions = []
        if not qualities:
            return regions

        current_region = {"start": 0, "high_quality": qualities[0] >= self.min_quality, "length": 1}

        for i in range(1, len(qualities)):
            is_high_quality = qualities[i] >= self.min_quality

            if is_high_quality == current_region["high_quality"]:
                current_region["length"] += 1
            else:
                current_region["end"] = current_region["start"] + current_region["length"] - 1
                regions.append(current_region)

                current_region = {"start": i, "high_quality": is_high_quality, "length": 1}

        # Add the last region
        current_region["end"] = current_region["start"] + current_region["length"] - 1
        regions.append(current_region)

        return regions
