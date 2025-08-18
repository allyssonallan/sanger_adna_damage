"""
Damage calculator component for ancient DNA analysis.
"""

import logging
from typing import Dict
from .sequence_aligner import SequenceAligner
from .damage_types import DamageAnalysisResult

logger = logging.getLogger(__name__)


class DamageCalculator:
    """Calculates damage statistics from aligned sequences."""
    
    def __init__(self, terminal_length: int = 10):
        """
        Initialize damage calculator.
        
        Args:
            terminal_length: Length of terminal regions to analyze for damage
        """
        self.terminal_length = terminal_length
        self.aligner = SequenceAligner()
        
    def calculate_damage_statistics(
        self, reference: str, query: str
    ) -> DamageAnalysisResult:
        """
        Calculate damage statistics from sequences.

        Args:
            reference: Reference sequence string
            query: Query sequence string

        Returns:
            Dictionary with damage statistics
        """
        # Align sequences
        alignment = self.aligner.align_sequences(reference, query)
        ref_clean, query_clean = self.aligner.clean_alignment(alignment[0], alignment[1])

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
        
    def calculate_positional_damage(
        self, reference: str, query: str
    ) -> Dict[str, list[float]]:
        """
        Calculate damage rates by position across the sequence.

        Args:
            reference: Reference sequence string
            query: Query sequence string

        Returns:
            Dictionary with positional damage rates
        """
        # Align sequences
        alignment = self.aligner.align_sequences(reference, query)
        ref_clean, query_clean = self.aligner.clean_alignment(alignment[0], alignment[1])

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
