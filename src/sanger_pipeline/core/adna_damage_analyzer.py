"""
Ancient DNA damage analysis module for detecting C->T deamination patterns.
"""

import logging
from pathlib import Path
from typing import Dict, List, Tuple
import numpy as np
import matplotlib.pyplot as plt
from Bio import SeqIO
from Bio.Align import PairwiseAligner

from ..utils.helpers import validate_file_exists


logger = logging.getLogger(__name__)


class ADNADamageAnalyzer:
    """
    Analyzer for ancient DNA damage patterns, specifically C->T deamination.
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

    def analyze_sequence_damage(
        self, sequence_file: Path, reference_file: Path
    ) -> Dict[str, float]:
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

        # Align sequences
        alignment = self._align_sequences(str(ref_seq.seq), str(query_seq.seq))

        # Calculate damage statistics
        damage_stats = self._calculate_damage_statistics(alignment)

        logger.info(f"Analyzed damage for {sequence_file.name}: "
                   f"5' damage = {damage_stats['damage_5_prime']:.3f}, "
                   f"3' damage = {damage_stats['damage_3_prime']:.3f}")

        return damage_stats

    def _align_sequences(self, reference: str, query: str) -> Tuple[str, str]:
        """
        Align two sequences using global alignment.

        Args:
            reference: Reference sequence string
            query: Query sequence string

        Returns:
            Tuple of (aligned_reference, aligned_query)
        """
        # Use Bio.Align.PairwiseAligner (modern replacement for pairwise2)
        aligner = PairwiseAligner()
        aligner.match_score = 2
        aligner.mismatch_score = -1
        aligner.open_gap_score = -2
        aligner.extend_gap_score = -0.5
        aligner.mode = 'global'
        
        alignments = aligner.align(reference, query)
        
        if not alignments:
            raise ValueError("No alignment found between sequences")
        
        # Take the best alignment (first one)
        best_alignment = alignments[0]
        aligned_ref = str(best_alignment[0])
        aligned_query = str(best_alignment[1])

        return aligned_ref, aligned_query

    def _calculate_damage_statistics(
        self, alignment: Tuple[str, str]
    ) -> Dict[str, float]:
        """
        Calculate damage statistics from alignment.

        Args:
            alignment: Tuple of (aligned_reference, aligned_query)

        Returns:
            Dictionary with damage statistics
        """
        aligned_ref, aligned_query = alignment

        # Remove gaps and get valid positions
        ref_clean, query_clean = self._clean_alignment(aligned_ref, aligned_query)

        # Calculate terminal damage
        damage_5_prime = self._calculate_terminal_damage(
            ref_clean, query_clean, "5_prime"
        )
        damage_3_prime = self._calculate_terminal_damage(
            ref_clean, query_clean, "3_prime"
        )

        # Calculate overall statistics
        total_bases = len(ref_clean)
        valid_bases = sum(1 for r, q in zip(ref_clean, query_clean) 
                         if r.upper() in "ACGT" and q.upper() in "ACGT")
        n_content = query_clean.upper().count("N")
        ambiguous_content = sum(1 for base in query_clean.upper() 
                               if base not in "ACGTN")
        
        total_ct_transitions = self._count_transitions(ref_clean, query_clean, "C", "T")
        total_ga_transitions = self._count_transitions(ref_clean, query_clean, "G", "A")

        # Calculate damage rate based on valid bases only
        overall_damage_rate = 0.0
        if valid_bases > 0:
            overall_damage_rate = (total_ct_transitions + total_ga_transitions) / valid_bases

        return {
            "damage_5_prime": damage_5_prime,
            "damage_3_prime": damage_3_prime,
            "total_bases": total_bases,
            "valid_bases": valid_bases,
            "n_content": n_content,
            "ambiguous_content": ambiguous_content,
            "sequence_quality": {
                "n_percentage": (n_content / total_bases * 100) if total_bases > 0 else 0,
                "valid_percentage": (valid_bases / total_bases * 100) if total_bases > 0 else 0
            },
            "total_ct_transitions": total_ct_transitions,
            "total_ga_transitions": total_ga_transitions,
            "overall_damage_rate": overall_damage_rate
        }

    def _clean_alignment(self, aligned_ref: str, aligned_query: str) -> Tuple[str, str]:
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

    def _calculate_terminal_damage(
        self, ref_seq: str, query_seq: str, terminal: str
    ) -> float:
        """
        Calculate damage rate at sequence terminals.

        Args:
            ref_seq: Reference sequence
            query_seq: Query sequence
            terminal: Either '5_prime' or '3_prime'

        Returns:
            Damage rate as float
        """
        seq_length = len(ref_seq)
        
        if seq_length < self.terminal_length:
            logger.warning(f"Sequence too short ({seq_length}) for terminal analysis")
            return 0.0

        if terminal == "5_prime":
            ref_region = ref_seq[:self.terminal_length]
            query_region = query_seq[:self.terminal_length]
        else:  # 3_prime
            ref_region = ref_seq[-self.terminal_length:]
            query_region = query_seq[-self.terminal_length:]

        # Count C->T and G->A transitions (excluding N bases)
        ct_transitions = self._count_transitions(ref_region, query_region, "C", "T")
        ga_transitions = self._count_transitions(ref_region, query_region, "G", "A")

        # Count valid positions where both ref and query have ACGT bases
        valid_positions = sum(1 for r, q in zip(ref_region, query_region) 
                             if r.upper() in "ACGT" and q.upper() in "ACGT")

        if valid_positions == 0:
            logger.warning(f"No valid bases in {terminal} terminal region for damage analysis")
            return 0.0

        damage_rate = (ct_transitions + ga_transitions) / valid_positions
        
        # Log quality metrics
        n_count = query_region.upper().count("N")
        logger.debug(f"{terminal} terminal: {valid_positions}/{len(ref_region)} valid positions, "
                    f"{n_count} N bases, damage rate: {damage_rate:.4f}")
        
        return damage_rate

    def _count_transitions(
        self, ref_seq: str, query_seq: str, from_base: str, to_base: str
    ) -> int:
        """
        Count specific base transitions, excluding N bases and ambiguous calls.

        Args:
            ref_seq: Reference sequence
            query_seq: Query sequence
            from_base: Base in reference
            to_base: Base in query

        Returns:
            Number of transitions
        """
        count = 0
        for ref_base, query_base in zip(ref_seq, query_seq):
            # Skip positions with N or ambiguous bases
            if (ref_base.upper() not in "ACGT" or 
                query_base.upper() not in "ACGT"):
                continue
            
            if ref_base.upper() == from_base and query_base.upper() == to_base:
                count += 1
        return count

    def bootstrap_damage_analysis(
        self, sequence_files: List[Path], reference_file: Path, 
        iterations: int = 10000
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

        # Perform bootstrap resampling
        bootstrap_5_prime = []
        bootstrap_3_prime = []

        for _ in range(iterations):
            # Resample with replacement
            sample_indices = np.random.choice(len(all_damage_data), 
                                            size=len(all_damage_data), 
                                            replace=True)
            
            sample_5_prime = [all_damage_data[i]["damage_5_prime"] for i in sample_indices]
            sample_3_prime = [all_damage_data[i]["damage_3_prime"] for i in sample_indices]
            
            bootstrap_5_prime.append(np.mean(sample_5_prime))
            bootstrap_3_prime.append(np.mean(sample_3_prime))

        # Calculate p-values (probability of observing damage <= 0.01)
        p_value_5_prime = np.mean(np.array(bootstrap_5_prime) <= 0.01)
        p_value_3_prime = np.mean(np.array(bootstrap_3_prime) <= 0.01)

        observed_5_prime = np.mean([data["damage_5_prime"] for data in all_damage_data])
        observed_3_prime = np.mean([data["damage_3_prime"] for data in all_damage_data])

        bootstrap_results = {
            "observed_damage_5_prime": observed_5_prime,
            "observed_damage_3_prime": observed_3_prime,
            "p_value_5_prime": p_value_5_prime,
            "p_value_3_prime": p_value_3_prime,
            "bootstrap_mean_5_prime": np.mean(bootstrap_5_prime),
            "bootstrap_mean_3_prime": np.mean(bootstrap_3_prime),
            "bootstrap_std_5_prime": np.std(bootstrap_5_prime),
            "bootstrap_std_3_prime": np.std(bootstrap_3_prime)
        }

        logger.info(f"Bootstrap analysis complete. "
                   f"5' damage: {observed_5_prime:.3f} (p={p_value_5_prime:.4f}), "
                   f"3' damage: {observed_3_prime:.3f} (p={p_value_3_prime:.4f})")

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
        
        damage_indicated_5 = (bootstrap_results["observed_damage_5_prime"] > self.min_damage_threshold and
                              bootstrap_results["p_value_5_prime"] < p_threshold)
        
        damage_indicated_3 = (bootstrap_results["observed_damage_3_prime"] > self.min_damage_threshold and
                              bootstrap_results["p_value_3_prime"] < p_threshold)

        overall_damage_indicated = damage_indicated_5 and damage_indicated_3

        if overall_damage_indicated:
            interpretation = ("Deamination patterns at both termini are consistent with aDNA damage signatures. "
                            "These patterns are indicative of ancient DNA, but additional validation through "
                            "independent extractions and contamination controls is recommended.")
            status = "DAMAGE_INDICATED"
        elif damage_indicated_5 or damage_indicated_3:
            interpretation = ("Modest deamination patterns detected. These may be indicative of some aDNA damage, "
                            "but the signal is not strong at both termini. Further validation recommended.")
            status = "PARTIAL_DAMAGE_SIGNATURE"
        else:
            interpretation = ("No significant terminal deamination excess detected. Patterns are consistent with "
                            "modern DNA or heavily degraded ancient DNA with low damage preservation.")
            status = "NO_DAMAGE_SIGNATURE"

        return {
            "status": status,
            "interpretation": interpretation,
            "5_prime_damage_indicated": bool(damage_indicated_5),
            "3_prime_damage_indicated": bool(damage_indicated_3)
        }

    def generate_damage_plots(
        self, sequence_files: List[Path], reference_file: Path, 
        output_dir: Path
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
                damage_profile = self._calculate_positional_damage(seq_file, reference_file)
                all_positional_damage.append({
                    "sample": seq_file.stem,
                    "profile": damage_profile
                })
            except Exception as e:
                logger.warning(f"Failed to generate damage profile for {seq_file}: {e}")
                continue

        if not all_positional_damage:
            logger.warning("No valid damage profiles generated")
            return

        # Create "smile plot" showing damage by position
        self._create_smile_plot(all_positional_damage, output_dir / "damage_smile_plot.png")
        
        # Create summary damage plot
        self._create_summary_damage_plot(all_positional_damage, output_dir / "damage_summary.png")

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

        # Align sequences
        alignment = self._align_sequences(str(ref_seq.seq), str(query_seq.seq))
        ref_clean, query_clean = self._clean_alignment(alignment[0], alignment[1])

        # Calculate damage rates for each position from terminals
        damage_5_prime = []
        damage_3_prime = []

        for i in range(min(self.terminal_length, len(ref_clean))):
            # 5' end
            if i < len(ref_clean):
                is_damage_5 = (ref_clean[i].upper() == "C" and query_clean[i].upper() == "T") or \
                             (ref_clean[i].upper() == "G" and query_clean[i].upper() == "A")
                damage_5_prime.append(1.0 if is_damage_5 else 0.0)

            # 3' end
            pos_3 = len(ref_clean) - 1 - i
            if pos_3 >= 0:
                is_damage_3 = (ref_clean[pos_3].upper() == "C" and query_clean[pos_3].upper() == "T") or \
                             (ref_clean[pos_3].upper() == "G" and query_clean[pos_3].upper() == "A")
                damage_3_prime.append(1.0 if is_damage_3 else 0.0)

        return {
            "5_prime": damage_5_prime,
            "3_prime": damage_3_prime
        }

    def _create_smile_plot(self, damage_data: List[Dict], output_file: Path) -> None:
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
                for i, rate in enumerate(profile["5_prime"][:self.terminal_length]):
                    avg_damage_5[i] += rate
                for i, rate in enumerate(profile["3_prime"][:self.terminal_length]):
                    avg_damage_3[i] += rate

        if valid_samples > 0:
            avg_damage_5 = [rate / valid_samples for rate in avg_damage_5]
            avg_damage_3 = [rate / valid_samples for rate in avg_damage_3]

        # Plot 5' damage
        ax.plot(positions_5, avg_damage_5, 'o-', color='red', 
               label="5' C→T + G→A", linewidth=2, markersize=6)
        
        # Plot 3' damage (reversed positions)
        ax.plot(positions_3, avg_damage_3, 'o-', color='blue', 
               label="3' C→T + G→A", linewidth=2, markersize=6)

        ax.set_xlabel("Distance from terminus (bp)")
        ax.set_ylabel("Misincorporation frequency")
        ax.set_title("aDNA Damage Pattern Analysis (Smile Plot)")
        ax.legend()
        ax.grid(True, alpha=0.3)
        ax.set_ylim(0, max(max(avg_damage_5 + avg_damage_3, default=0.1) * 1.1, 0.05))

        plt.tight_layout()
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        plt.close()

    def _create_summary_damage_plot(self, damage_data: List[Dict], output_file: Path) -> None:
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
        ax1.bar(range(len(samples)), damage_5, color='red', alpha=0.7)
        ax1.axhline(y=self.min_damage_threshold, color='black', linestyle='--', 
                   label=f'Threshold ({self.min_damage_threshold})')
        ax1.set_title("5' Terminal Damage")
        ax1.set_ylabel("Damage Rate")
        ax1.set_xticks(range(len(samples)))
        ax1.set_xticklabels(samples, rotation=45, ha='right')
        ax1.legend()

        # 3' damage plot
        ax2.bar(range(len(samples)), damage_3, color='blue', alpha=0.7)
        ax2.axhline(y=self.min_damage_threshold, color='black', linestyle='--', 
                   label=f'Threshold ({self.min_damage_threshold})')
        ax2.set_title("3' Terminal Damage")
        ax2.set_ylabel("Damage Rate")
        ax2.set_xticks(range(len(samples)))
        ax2.set_xticklabels(samples, rotation=45, ha='right')
        ax2.legend()

        plt.tight_layout()
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        plt.close()
