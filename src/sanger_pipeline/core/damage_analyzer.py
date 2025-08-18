"""
Ancient DNA damage analysis components for Enhanced AB1 Converter.

This module provides specialized functionality for analyzing and handling
ancient DNA damage patterns in sequences.
"""

import logging
from typing import Dict

logger = logging.getLogger(__name__)


class DamageAnalyzer:
    """Analyzes ancient DNA damage patterns in sequences."""

    def __init__(self, extremity_analysis_length: int = 30):
        """
        Initialize damage analyzer.

        Args:
            extremity_analysis_length: Length of sequence extremities to analyze
        """
        self.extremity_analysis_length = extremity_analysis_length

    def analyze_extremity_n_abundance(
        self, sequence: str, stage: str = "unknown"
    ) -> Dict:
        """
        Analyze N abundance at sequence extremities (aDNA damage indicator).

        Ancient DNA often shows C→T and G→A transitions at 5' and 3' ends
        respectively, which can appear as N calls in poor quality regions.

        Args:
            sequence: DNA sequence to analyze
            stage: Processing stage for logging

        Returns:
            Dictionary with N abundance analysis
        """
        if not sequence:
            return {
                "five_prime_n_percentage": 0.0,
                "three_prime_n_percentage": 0.0,
                "total_n_percentage": 0.0,
                "damage_pattern": "no_sequence",
                "analysis_length": 0,
                "stage": stage,
            }

        seq_length = len(sequence)
        analysis_length = min(self.extremity_analysis_length, seq_length // 2)

        if analysis_length <= 0:
            return {
                "five_prime_n_percentage": 0.0,
                "three_prime_n_percentage": 0.0,
                "total_n_percentage": 0.0,
                "damage_pattern": "sequence_too_short",
                "analysis_length": 0,
                "stage": stage,
            }

        # Analyze 5' end
        five_prime_region = sequence[:analysis_length].upper()
        five_prime_n_count = five_prime_region.count("N")
        five_prime_n_percentage = (five_prime_n_count / analysis_length) * 100

        # Analyze 3' end
        three_prime_region = sequence[-analysis_length:].upper()
        three_prime_n_count = three_prime_region.count("N")
        three_prime_n_percentage = (three_prime_n_count / analysis_length) * 100

        # Total N percentage in whole sequence
        total_n_count = sequence.upper().count("N")
        total_n_percentage = (total_n_count / seq_length) * 100

        # Classify damage pattern
        damage_pattern = self._classify_damage_pattern(
            five_prime_n_percentage, three_prime_n_percentage, total_n_percentage
        )

        return {
            "five_prime_n_percentage": round(five_prime_n_percentage, 2),
            "three_prime_n_percentage": round(three_prime_n_percentage, 2),
            "total_n_percentage": round(total_n_percentage, 2),
            "damage_pattern": damage_pattern,
            "analysis_length": analysis_length,
            "stage": stage,
            "five_prime_region": five_prime_region,
            "three_prime_region": three_prime_region,
            "sequence_length": seq_length,
        }

    def _classify_damage_pattern(
        self, five_prime_n: float, three_prime_n: float, total_n: float
    ) -> str:
        """
        Classify damage pattern based on N abundance.

        Args:
            five_prime_n: N percentage at 5' end
            three_prime_n: N percentage at 3' end
            total_n: Total N percentage

        Returns:
            Damage pattern classification
        """
        # Define thresholds for damage classification
        moderate_threshold = 15.0
        high_threshold = 30.0
        extreme_threshold = 50.0

        # Check for extreme damage patterns
        if five_prime_n > extreme_threshold or three_prime_n > extreme_threshold:
            return "extreme_damage"

        # Check for high damage at both ends (typical aDNA)
        if five_prime_n > high_threshold and three_prime_n > high_threshold:
            return "high_bilateral_damage"

        # Check for high damage at one end
        if five_prime_n > high_threshold or three_prime_n > high_threshold:
            if five_prime_n > three_prime_n:
                return "high_five_prime_damage"
            else:
                return "high_three_prime_damage"

        # Check for moderate damage patterns
        if five_prime_n > moderate_threshold or three_prime_n > moderate_threshold:
            return "moderate_damage"

        # Check for asymmetric damage (more common in aDNA)
        if abs(five_prime_n - three_prime_n) > 10.0:
            return "asymmetric_damage"

        # Low or no damage
        if total_n < 5.0:
            return "minimal_damage"
        else:
            return "distributed_damage"

    def adjust_parameters_for_damage(self, n_analysis: Dict) -> Dict:
        """
        Suggest parameter adjustments based on damage analysis.

        Args:
            n_analysis: N abundance analysis results

        Returns:
            Dictionary with suggested parameter adjustments
        """
        damage_pattern = n_analysis.get("damage_pattern", "unknown")
        # Note: Individual N percentages extracted but used in pattern-based logic

        adjustments = {
            "quality_adjustment": 0,
            "trim_adjustment": 0,
            "min_length_adjustment": 0,
            "reasoning": "",
            "severity": "none",
        }

        if damage_pattern == "extreme_damage":
            adjustments.update(
                {
                    "quality_adjustment": -5,  # Lower quality threshold
                    "trim_adjustment": 10,  # More aggressive trimming
                    "min_length_adjustment": -20,  # Accept shorter sequences
                    "reasoning": "Extreme damage detected - using relaxed parameters",
                    "severity": "extreme",
                }
            )

        elif damage_pattern in [
            "high_bilateral_damage",
            "high_five_prime_damage",
            "high_three_prime_damage",
        ]:
            adjustments.update(
                {
                    "quality_adjustment": -3,
                    "trim_adjustment": 5,
                    "min_length_adjustment": -10,
                    "reasoning": "High damage at sequence extremities - moderately relaxed parameters",
                    "severity": "high",
                }
            )

        elif damage_pattern == "moderate_damage":
            adjustments.update(
                {
                    "quality_adjustment": -2,
                    "trim_adjustment": 3,
                    "min_length_adjustment": -5,
                    "reasoning": "Moderate damage detected - slightly relaxed parameters",
                    "severity": "moderate",
                }
            )

        elif damage_pattern == "asymmetric_damage":
            adjustments.update(
                {
                    "quality_adjustment": -1,
                    "trim_adjustment": 2,
                    "min_length_adjustment": 0,
                    "reasoning": "Asymmetric damage pattern - minor adjustments",
                    "severity": "low",
                }
            )

        else:  # minimal_damage, distributed_damage, no damage
            adjustments.update(
                {
                    "reasoning": "No significant damage detected - using standard parameters",
                    "severity": "none",
                }
            )

        return adjustments

    def create_damage_summary(self, processing_stats: Dict) -> Dict:
        """
        Create comprehensive damage analysis summary.

        Args:
            processing_stats: Processing statistics from conversion

        Returns:
            Damage analysis summary
        """
        stages = ["original", "after_primer_removal", "after_quality_filter", "final"]

        summary = {
            "damage_progression": {},
            "overall_assessment": "",
            "recommendations": [],
            "damage_indicators": {},
        }

        # Analyze damage progression through stages
        for stage in stages:
            stage_key = f"{stage}_n_analysis"
            if stage_key in processing_stats:
                n_analysis = processing_stats[stage_key]
                summary["damage_progression"][stage] = {
                    "five_prime_damage": n_analysis.get("five_prime_n_percentage", 0.0),
                    "three_prime_damage": n_analysis.get(
                        "three_prime_n_percentage", 0.0
                    ),
                    "total_damage": n_analysis.get("total_n_percentage", 0.0),
                    "pattern": n_analysis.get("damage_pattern", "unknown"),
                }

        # Determine overall assessment
        if processing_stats.get("final_n_analysis", {}).get("damage_pattern", "") in [
            "extreme_damage",
            "high_bilateral_damage",
        ]:
            summary["overall_assessment"] = "high_damage"
            summary["recommendations"] = [
                "Consider using more relaxed quality parameters",
                "Validate results with independent methods",
                "Check for contamination or degradation",
            ]
        elif processing_stats.get("final_n_analysis", {}).get("damage_pattern", "") in [
            "moderate_damage",
            "high_five_prime_damage",
            "high_three_prime_damage",
        ]:
            summary["overall_assessment"] = "moderate_damage"
            summary["recommendations"] = [
                "Results appear reasonable for ancient DNA",
                "Consider damage correction in downstream analysis",
            ]
        else:
            summary["overall_assessment"] = "low_damage"
            summary["recommendations"] = [
                "High quality sequence with minimal damage",
                "Suitable for standard analysis workflows",
            ]

        # Compile damage indicators
        final_analysis = processing_stats.get("final_n_analysis", {})
        summary["damage_indicators"] = {
            "ancient_dna_compatible": final_analysis.get("damage_pattern", "")
            not in ["extreme_damage"],
            "quality_concerns": final_analysis.get("total_n_percentage", 0.0) > 20.0,
            "asymmetric_pattern": "asymmetric"
            in final_analysis.get("damage_pattern", ""),
            "bilateral_damage": "bilateral" in final_analysis.get("damage_pattern", ""),
        }

        return summary
