"""
Type definitions for damage analysis components.
"""

from typing import TypedDict


class SequenceQuality(TypedDict):
    """Quality metrics for a sequence."""

    n_percentage: float
    valid_percentage: float


class DamageAnalysisResult(TypedDict):
    """Results from damage analysis."""

    damage_5_prime: float
    damage_3_prime: float
    total_bases: int
    valid_bases: int
    n_content: int
    ambiguous_content: int
    sequence_quality: SequenceQuality
    total_ct_transitions: int
    total_ga_transitions: int
    overall_damage_rate: float
