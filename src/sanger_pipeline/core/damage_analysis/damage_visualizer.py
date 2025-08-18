"""
Visualization component for damage analysis plots.
"""

import logging
from pathlib import Path
from typing import Dict, List
import numpy as np
import matplotlib.pyplot as plt

logger = logging.getLogger(__name__)


class DamageVisualizer:
    """Generates damage pattern visualization plots."""

    def __init__(self, terminal_length: int = 10, min_damage_threshold: float = 0.02):
        """
        Initialize damage visualizer.

        Args:
            terminal_length: Length of terminal regions to analyze
            min_damage_threshold: Minimum damage threshold for plots
        """
        self.terminal_length = terminal_length
        self.min_damage_threshold = min_damage_threshold

    def create_smile_plot(self, damage_data: List[Dict], output_file: Path) -> None:
        """
        Create the characteristic "smile plot" for aDNA damage.

        Args:
            damage_data: List of damage profiles
            output_file: Output file path
        """
        fig, ax = plt.subplots(figsize=(12, 6))

        # Aggregate damage rates by position
        positions_5 = list(range(1, self.terminal_length + 1))
        positions_3 = list(range(1, self.terminal_length + 1))

        avg_damage_5 = [0] * self.terminal_length
        avg_damage_3 = [0] * self.terminal_length

        valid_samples = 0

        for sample_data in damage_data:
            profile = sample_data["profile"]
            if "5_prime" in profile and "3_prime" in profile:
                valid_samples += 1
                for i, rate in enumerate(profile["5_prime"][: self.terminal_length]):
                    avg_damage_5[i] += rate
                for i, rate in enumerate(profile["3_prime"][: self.terminal_length]):
                    avg_damage_3[i] += rate

        if valid_samples > 0:
            avg_damage_5 = [rate / valid_samples for rate in avg_damage_5]
            avg_damage_3 = [rate / valid_samples for rate in avg_damage_3]

        # Plot 5' damage
        ax.plot(
            positions_5,
            avg_damage_5,
            "o-",
            color="red",
            label="5' C→T + G→A",
            linewidth=2,
            markersize=6,
        )

        # Plot 3' damage (reversed positions)
        ax.plot(
            positions_3,
            avg_damage_3,
            "o-",
            color="blue",
            label="3' C→T + G→A",
            linewidth=2,
            markersize=6,
        )

        ax.set_xlabel("Distance from terminus (bp)")
        ax.set_ylabel("Misincorporation frequency")
        ax.set_title("aDNA Damage Pattern Analysis (Smile Plot)")
        ax.legend()
        ax.grid(True, alpha=0.3)
        ax.set_ylim(0, max(max(avg_damage_5 + avg_damage_3, default=0.1) * 1.1, 0.05))

        plt.tight_layout()
        plt.savefig(output_file, dpi=300, bbox_inches="tight")
        plt.close()

    def create_summary_damage_plot(
        self, damage_data: List[Dict], output_file: Path
    ) -> None:
        """
        Create summary damage plot for all samples.

        Args:
            damage_data: List of damage profiles
            output_file: Output file path
        """
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))

        samples = [data["sample"] for data in damage_data]
        damage_5 = []
        damage_3 = []

        for sample_data in damage_data:
            profile = sample_data["profile"]
            # Calculate average damage for each terminal
            avg_5 = np.mean(profile.get("5_prime", [0]))
            avg_3 = np.mean(profile.get("3_prime", [0]))
            damage_5.append(avg_5)
            damage_3.append(avg_3)

        # 5' damage plot
        ax1.bar(range(len(samples)), damage_5, color="red", alpha=0.7)
        ax1.axhline(
            y=self.min_damage_threshold,
            color="black",
            linestyle="--",
            label=f"Threshold ({self.min_damage_threshold})",
        )
        ax1.set_title("5' Terminal Damage")
        ax1.set_ylabel("Damage Rate")
        ax1.set_xticks(range(len(samples)))
        ax1.set_xticklabels(samples, rotation=45, ha="right")
        ax1.legend()

        # 3' damage plot
        ax2.bar(range(len(samples)), damage_3, color="blue", alpha=0.7)
        ax2.axhline(
            y=self.min_damage_threshold,
            color="black",
            linestyle="--",
            label=f"Threshold ({self.min_damage_threshold})",
        )
        ax2.set_title("3' Terminal Damage")
        ax2.set_ylabel("Damage Rate")
        ax2.set_xticks(range(len(samples)))
        ax2.set_xticklabels(samples, rotation=45, ha="right")
        ax2.legend()

        plt.tight_layout()
        plt.savefig(output_file, dpi=300, bbox_inches="tight")
        plt.close()
