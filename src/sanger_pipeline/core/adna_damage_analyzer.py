"""
Refactored ancient DNA damage analysis module using modular components.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Dict, List
from Bio import SeqIO

from ..utils.helpers import validate_file_exists
from .damage_analysis import (
    DamageAnalysisResult,
    DamageCalculator,
    StatisticalAnalyzer,
    DamageVisualizer,
)

logger = logging.getLogger(__name__)


class ADNADamageAnalyzer:
    """
    Refactored analyzer for ancient DNA damage patterns using modular components.
    """

    def __init__(self, terminal_length: int = 10, min_damage_threshold: float = 0.02):
        """
        Initialize aDNA damage analyzer.

        Args:
            terminal_length: Length of terminal regions to analyze for damage
            min_damage_threshold: Minimum damage rate to consider significant
        """
        self.terminal_length = terminal_length
        self.min_damage_threshold = min_damage_threshold

        # Initialize components
        self.damage_calculator = DamageCalculator(terminal_length)
        self.statistical_analyzer = StatisticalAnalyzer(min_damage_threshold)
        self.visualizer = DamageVisualizer(terminal_length, min_damage_threshold)

    def analyze_sequence_damage(
        self, sequence_file: Path, reference_file: Path
    ) -> DamageAnalysisResult:
        """
        Analyze damage patterns in a sequence against reference.

        Args:
            sequence_file: Path to query sequence FASTA file
            reference_file: Path to reference sequence FASTA file

        Returns:
            Dictionary with damage statistics
        """
        validate_file_exists(sequence_file, "sequence file")
        validate_file_exists(reference_file, "reference file")

        # Read sequences
        query_seq = SeqIO.read(sequence_file, "fasta")
        ref_seq = SeqIO.read(reference_file, "fasta")

        # Calculate damage statistics using modular component
        damage_stats = self.damage_calculator.calculate_damage_statistics(
            str(ref_seq.seq), str(query_seq.seq)
        )

        logger.info(
            f"Analyzed damage for {sequence_file.name}: "
            f"5' damage = {damage_stats['damage_5_prime']:.3f}, "
            f"3' damage = {damage_stats['damage_3_prime']:.3f}"
        )

        return damage_stats

    def bootstrap_damage_analysis(
        self, sequence_files: List[Path], reference_file: Path, iterations: int = 10000
    ) -> Dict:
        """
        Perform bootstrap analysis for damage assessment.

        Args:
            sequence_files: List of sequence files to analyze
            reference_file: Reference sequence file
            iterations: Number of bootstrap iterations

        Returns:
            Bootstrap statistics
        """
        logger.info(f"Starting bootstrap analysis with {iterations} iterations")

        # Collect damage data from all sequences
        all_damage_data = []
        for seq_file in sequence_files:
            try:
                damage_stats = self.analyze_sequence_damage(seq_file, reference_file)
                all_damage_data.append(damage_stats)
            except Exception as e:
                logger.warning(f"Failed to analyze {seq_file}: {e}")
                continue

        if not all_damage_data:
            raise ValueError("No valid sequences for bootstrap analysis")

        # Use statistical analyzer component
        return self.statistical_analyzer.bootstrap_damage_analysis(
            all_damage_data, iterations
        )

    def assess_damage_indicators(self, bootstrap_results: Dict) -> Dict:
        """
        Assess aDNA damage indicators based on deamination patterns.

        Args:
            bootstrap_results: Results from bootstrap analysis

        Returns:
            Damage pattern assessment indicating potential aDNA signatures
        """
        return self.statistical_analyzer.assess_damage_indicators(bootstrap_results)

    def assess_authenticity(self, bootstrap_results: Dict) -> Dict:
        """
        Assess aDNA authenticity (backward compatibility alias).

        Args:
            bootstrap_results: Results from bootstrap analysis

        Returns:
            Damage pattern assessment indicating potential aDNA signatures
        """
        return self.assess_damage_indicators(bootstrap_results)

    def generate_damage_plots(
        self, sequence_files: List[Path], reference_file: Path, output_dir: Path
    ) -> None:
        """
        Generate damage pattern visualization plots.

        Args:
            sequence_files: List of sequence files
            reference_file: Reference sequence file
            output_dir: Directory to save plots
        """
        output_dir.mkdir(parents=True, exist_ok=True)

        # Analyze all sequences for positional damage
        all_positional_damage = []

        for seq_file in sequence_files:
            try:
                damage_profile = self._calculate_positional_damage(
                    seq_file, reference_file
                )
                all_positional_damage.append(
                    {"sample": seq_file.stem, "profile": damage_profile}
                )
            except Exception as e:
                logger.warning(f"Failed to generate damage profile for {seq_file}: {e}")
                continue

        if not all_positional_damage:
            logger.warning("No valid damage profiles generated")
            return

        # Use visualizer component to create plots
        self.visualizer.create_smile_plot(
            all_positional_damage, output_dir / "damage_smile_plot.png"
        )

        self.visualizer.create_summary_damage_plot(
            all_positional_damage, output_dir / "damage_summary.png"
        )

        logger.info(f"Damage plots saved to {output_dir}")

    def _calculate_positional_damage(
        self, sequence_file: Path, reference_file: Path
    ) -> Dict[str, List[float]]:
        """
        Calculate damage rates by position across the sequence.

        Args:
            sequence_file: Sequence file to analyze
            reference_file: Reference sequence file

        Returns:
            Dictionary with positional damage rates
        """
        query_seq = SeqIO.read(sequence_file, "fasta")
        ref_seq = SeqIO.read(reference_file, "fasta")

        # Use damage calculator component
        return self.damage_calculator.calculate_positional_damage(
            str(ref_seq.seq), str(query_seq.seq)
        )
