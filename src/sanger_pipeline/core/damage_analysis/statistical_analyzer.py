"""
Statistical analysis component for damage assessment.
"""

import logging
from typing import Dict, List
import numpy as np

logger = logging.getLogger(__name__)


class StatisticalAnalyzer:
    """Handles bootstrap analysis and statistical assessment of damage patterns."""

    def __init__(self, min_damage_threshold: float = 0.02):
        """
        Initialize statistical analyzer.

        Args:
            min_damage_threshold: Minimum damage rate to consider significant
        """
        self.min_damage_threshold = min_damage_threshold

    def bootstrap_damage_analysis(
        self, damage_data: List[Dict], iterations: int = 10000
    ) -> Dict:
        """
        Perform bootstrap analysis for damage assessment.

        Args:
            damage_data: List of damage analysis results
            iterations: Number of bootstrap iterations

        Returns:
            Bootstrap statistics
        """
        logger.info(f"Starting bootstrap analysis with {iterations} iterations")

        if not damage_data:
            raise ValueError("No valid sequences for bootstrap analysis")

        # Perform bootstrap resampling
        bootstrap_5_prime = []
        bootstrap_3_prime = []

        for _ in range(iterations):
            # Resample with replacement
            sample_indices = np.random.choice(
                len(damage_data), size=len(damage_data), replace=True
            )

            sample_5_prime = [damage_data[i]["damage_5_prime"] for i in sample_indices]
            sample_3_prime = [damage_data[i]["damage_3_prime"] for i in sample_indices]

            bootstrap_5_prime.append(np.mean(sample_5_prime))
            bootstrap_3_prime.append(np.mean(sample_3_prime))

        # Calculate p-values (probability of observing damage <= 0.01)
        p_value_5_prime = np.mean(np.array(bootstrap_5_prime) <= 0.01)
        p_value_3_prime = np.mean(np.array(bootstrap_3_prime) <= 0.01)

        observed_5_prime = np.mean([data["damage_5_prime"] for data in damage_data])
        observed_3_prime = np.mean([data["damage_3_prime"] for data in damage_data])

        bootstrap_results = {
            "observed_damage_5_prime": observed_5_prime,
            "observed_damage_3_prime": observed_3_prime,
            "p_value_5_prime": p_value_5_prime,
            "p_value_3_prime": p_value_3_prime,
            "bootstrap_mean_5_prime": np.mean(bootstrap_5_prime),
            "bootstrap_mean_3_prime": np.mean(bootstrap_3_prime),
            "bootstrap_std_5_prime": np.std(bootstrap_5_prime),
            "bootstrap_std_3_prime": np.std(bootstrap_3_prime),
        }

        logger.info(
            f"Bootstrap analysis complete. "
            f"5' damage: {observed_5_prime:.3f} (p={p_value_5_prime:.4f}), "
            f"3' damage: {observed_3_prime:.3f} (p={p_value_3_prime:.4f})"
        )

        return bootstrap_results

    def assess_damage_indicators(self, bootstrap_results: Dict) -> Dict:
        """
        Assess aDNA damage indicators based on deamination patterns.

        Args:
            bootstrap_results: Results from bootstrap analysis

        Returns:
            Damage pattern assessment indicating potential aDNA signatures
        """
        p_threshold = 0.05

        damage_indicated_5 = (
            bootstrap_results["observed_damage_5_prime"] > self.min_damage_threshold
            and bootstrap_results["p_value_5_prime"] < p_threshold
        )

        damage_indicated_3 = (
            bootstrap_results["observed_damage_3_prime"] > self.min_damage_threshold
            and bootstrap_results["p_value_3_prime"] < p_threshold
        )

        overall_damage_indicated = damage_indicated_5 and damage_indicated_3

        if overall_damage_indicated:
            interpretation = (
                "Deamination patterns at both termini are consistent with aDNA damage signatures. "
                "These patterns are indicative of ancient DNA, but additional validation through "
                "independent extractions and contamination controls is recommended."
            )
            status = "DAMAGE_INDICATED"
        elif damage_indicated_5 or damage_indicated_3:
            interpretation = (
                "Modest deamination patterns detected. These may be indicative of some aDNA damage, "
                "but the signal is not strong at both termini. Further validation recommended."
            )
            status = "PARTIAL_DAMAGE_SIGNATURE"
        else:
            interpretation = (
                "No significant terminal deamination excess detected. Patterns are consistent with "
                "modern DNA or heavily degraded ancient DNA with low damage preservation."
            )
            status = "NO_DAMAGE_SIGNATURE"

        return {
            "status": status,
            "interpretation": interpretation,
            "5_prime_damage_indicated": bool(damage_indicated_5),
            "3_prime_damage_indicated": bool(damage_indicated_3),
        }
