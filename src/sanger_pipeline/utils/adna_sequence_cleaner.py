#!/usr/bin/env python3
"""
Ancient DNA Sequence Cleaner and Quality Control

This module provides tools for cleaning and quality control of ancient DNA sequences
before HSD conversion, specifically addressing common aDNA artifacts.

Author: Sanger aDNA Pipeline
"""

import logging
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import argparse

logger = logging.getLogger(__name__)


class aDNASequenceCleaner:
    """
    Sequence cleaner specifically designed for ancient DNA artifacts.
    """

    def __init__(self, min_length: int = 50, min_quality: float = 0.6):
        """
        Initialize the sequence cleaner.

        Args:
            min_length: Minimum sequence length to keep
            min_quality: Minimum proportion of valid nucleotides
        """
        self.min_length = min_length
        self.min_quality = min_quality

        # Valid DNA nucleotides
        self.valid_bases = {"A", "T", "G", "C"}

        # Common aDNA damage patterns
        self.damage_patterns = {
            "C->T": ("C", "T"),  # 5' C to T transitions
            "G->A": ("G", "A"),  # 3' G to A transitions
        }

        # Ambiguous nucleotide resolution priorities
        self.ambiguous_resolution = {
            "R": "G",  # Purine -> prefer G
            "Y": "C",  # Pyrimidine -> prefer C
            "S": "G",  # Strong -> prefer G
            "W": "A",  # Weak -> prefer A
            "K": "G",  # Keto -> prefer G
            "M": "A",  # Amino -> prefer A
            "B": "G",  # Not A -> prefer G
            "D": "G",  # Not C -> prefer G
            "H": "A",  # Not G -> prefer A
            "V": "G",  # Not T -> prefer G
            "N": "N",  # Unknown -> keep as N
        }

    def calculate_quality_score(self, sequence: str) -> float:
        """
        Calculate sequence quality based on valid nucleotides.

        Args:
            sequence: DNA sequence string

        Returns:
            Quality score (0.0 to 1.0)
        """
        if not sequence:
            return 0.0

        valid_count = sum(1 for base in sequence.upper() if base in self.valid_bases)
        return valid_count / len(sequence)

    def clean_sequence(self, sequence: str, aggressive_cleaning: bool = False) -> str:
        """
        Clean a single DNA sequence.

        Args:
            sequence: Input DNA sequence
            aggressive_cleaning: Apply more aggressive cleaning

        Returns:
            Cleaned sequence
        """
        cleaned = sequence.upper()

        # Step 1: Resolve ambiguous nucleotides
        cleaned = self._resolve_ambiguous_nucleotides(cleaned)

        # Step 2: Filter low quality regions
        if aggressive_cleaning:
            cleaned = self._filter_low_quality_windows(cleaned)

        # Step 3: Remove leading/trailing Ns
        cleaned = cleaned.strip("N")

        return cleaned

    def _resolve_ambiguous_nucleotides(self, sequence: str) -> str:
        """
        Resolve ambiguous nucleotides using priority rules.

        Args:
            sequence: Input sequence

        Returns:
            Sequence with resolved ambiguous nucleotides
        """
        resolved = []

        for base in sequence:
            if base in self.valid_bases:
                resolved.append(base)
            elif base in self.ambiguous_resolution:
                resolved.append(self.ambiguous_resolution[base])
            else:
                resolved.append("N")

        return "".join(resolved)

    def _filter_low_quality_windows(self, sequence: str, window_size: int = 15) -> str:
        """
        Filter low quality windows in the sequence.

        Args:
            sequence: Input sequence
            window_size: Window size for quality assessment

        Returns:
            Filtered sequence
        """
        if len(sequence) < window_size:
            return sequence

        filtered = list(sequence)

        for i in range(len(sequence) - window_size + 1):
            window = sequence[i : i + window_size]
            window_quality = self.calculate_quality_score(window)

            if window_quality < self.min_quality:
                # Mark this window as low quality
                for j in range(i, min(i + window_size, len(filtered))):
                    filtered[j] = "N"

        return "".join(filtered)

    def _remove_poly_n_regions(self, sequence: str, min_poly_n: int = 5) -> str:
        """
        Remove regions with excessive Ns.

        Args:
            sequence: Input sequence
            min_poly_n: Minimum consecutive Ns to remove

        Returns:
            Sequence with poly-N regions removed
        """
        import re

        # Replace regions with 5+ consecutive Ns with a single N
        pattern = f"N{{{min_poly_n},}}"
        return re.sub(pattern, "N", sequence)

    def clean_fasta_file(
        self, input_file: str, output_file: str, aggressive_cleaning: bool = False
    ) -> None:
        """
        Clean all sequences in a FASTA file.

        Args:
            input_file: Path to input FASTA file
            output_file: Path to output cleaned FASTA file
            aggressive_cleaning: Apply aggressive cleaning
        """
        logger.info(f"Cleaning FASTA file: {input_file}")
        logger.info(f"Minimum length: {self.min_length}")
        logger.info(f"Minimum quality: {self.min_quality}")
        logger.info(f"Aggressive cleaning: {aggressive_cleaning}")

        cleaned_records = []
        total_sequences = 0
        kept_sequences = 0

        for record in SeqIO.parse(input_file, "fasta"):
            total_sequences += 1
            original_seq = str(record.seq)

            # Calculate original quality
            original_quality = self.calculate_quality_score(original_seq)

            # Clean the sequence
            cleaned_seq = self.clean_sequence(original_seq, aggressive_cleaning)

            # Calculate cleaned quality
            cleaned_quality = self.calculate_quality_score(cleaned_seq)

            # Apply filters
            if (
                len(cleaned_seq) >= self.min_length
                and cleaned_quality >= self.min_quality
            ):

                # Create new record with cleaned sequence
                new_record = SeqRecord(
                    Seq(cleaned_seq),
                    id=record.id,
                    description=f"{record.description} | Quality: {cleaned_quality:.3f} | Length: {len(cleaned_seq)}",
                )

                cleaned_records.append(new_record)
                kept_sequences += 1

                logger.info(
                    f"Kept {record.id}: {len(original_seq)}bp -> {len(cleaned_seq)}bp, "
                    f"Quality: {original_quality:.3f} -> {cleaned_quality:.3f}"
                )
            else:
                logger.warning(
                    f"Filtered {record.id}: Length={len(cleaned_seq)}, "
                    f"Quality={cleaned_quality:.3f}"
                )

        # Write cleaned sequences
        with open(output_file, "w") as f:
            SeqIO.write(cleaned_records, f, "fasta")

        logger.info(
            f"Cleaning completed: {kept_sequences}/{total_sequences} sequences kept"
        )
        logger.info(f"Output written to: {output_file}")


def main():
    """Main function for command line usage."""
    parser = argparse.ArgumentParser(description="Ancient DNA Sequence Cleaner")
    parser.add_argument("-i", "--input", required=True, help="Input FASTA file")
    parser.add_argument(
        "-o", "--output", required=True, help="Output cleaned FASTA file"
    )
    parser.add_argument(
        "-l", "--min-length", type=int, default=50, help="Minimum sequence length"
    )
    parser.add_argument(
        "-q",
        "--min-quality",
        type=float,
        default=0.6,
        help="Minimum quality score (0.0-1.0)",
    )
    parser.add_argument(
        "--aggressive", action="store_true", help="Apply aggressive cleaning"
    )
    parser.add_argument("-v", "--verbose", action="store_true", help="Verbose logging")

    args = parser.parse_args()

    # Set up logging
    log_level = logging.DEBUG if args.verbose else logging.INFO
    logging.basicConfig(
        level=log_level, format="%(asctime)s - %(levelname)s - %(message)s"
    )

    # Create cleaner and process
    cleaner = aDNASequenceCleaner(
        min_length=args.min_length, min_quality=args.min_quality
    )

    cleaner.clean_fasta_file(
        input_file=args.input,
        output_file=args.output,
        aggressive_cleaning=args.aggressive,
    )


if __name__ == "__main__":
    main()
