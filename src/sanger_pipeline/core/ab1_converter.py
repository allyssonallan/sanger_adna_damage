"""
AB1 file converter with quality filtering and plotting capabilities.
"""

import logging
from pathlib import Path
from typing import Tuple, Optional
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import matplotlib.pyplot as plt

from ..utils.constants import DEFAULT_MIN_QUALITY, DEFAULT_MIN_SEQUENCE_LENGTH
from ..utils.helpers import validate_file_exists


logger = logging.getLogger(__name__)


class AB1Converter:
    """
    Converter for AB1 files to FASTA format with quality filtering.
    """

    def __init__(
        self,
        min_quality: int = DEFAULT_MIN_QUALITY,
        min_sequence_length: int = DEFAULT_MIN_SEQUENCE_LENGTH,
    ):
        """
        Initialize AB1 converter.

        Args:
            min_quality: Minimum Phred quality score for base calling
            min_sequence_length: Minimum sequence length after filtering
        """
        self.min_quality = min_quality
        self.min_sequence_length = min_sequence_length

    def convert_to_fasta(self, ab1_file: Path, output_file: Path) -> SeqRecord:
        """
        Convert AB1 file to FASTA format.

        Args:
            ab1_file: Path to input AB1 file
            output_file: Path to output FASTA file

        Returns:
            SeqRecord object of the converted sequence

        Raises:
            FileNotFoundError: If AB1 file doesn't exist
            ValueError: If AB1 file cannot be parsed
        """
        validate_file_exists(ab1_file, "AB1 file")

        try:
            record = SeqIO.read(ab1_file, "abi")
            logger.info(f"Read AB1 file: {ab1_file}")
        except Exception as e:
            raise ValueError(f"Failed to parse AB1 file {ab1_file}: {e}")

        # Set FASTA header to match output filename (with .fasta extension)
        record.id = output_file.name
        record.description = ""

        # Ensure output directory exists
        output_file.parent.mkdir(parents=True, exist_ok=True)

        # Write raw FASTA
        SeqIO.write(record, output_file, "fasta")
        logger.info(f"Wrote FASTA file: {output_file}")

        return record

    def validate_sequence_length(self, sequence: str) -> bool:
        """
        Check if sequence meets minimum length requirement after removing N's and gaps.

        Args:
            sequence: DNA sequence string

        Returns:
            True if sequence is long enough, False otherwise
        """
        # Count only valid nucleotides (not N's or gaps)
        valid_bases = sum(1 for base in sequence if base in "ATCG")
        is_valid = valid_bases >= self.min_sequence_length

        if not is_valid:
            logger.warning(
                f"Sequence too short: {valid_bases} valid bases (minimum: {self.min_sequence_length})"
            )

        return is_valid

    def filter_by_quality(
        self, record: SeqRecord, output_file: Path
    ) -> Optional[SeqRecord]:
        """
        Filter sequence by quality scores, replacing low-quality bases with 'N',
        and check minimum sequence length.

        Args:
            record: Input SeqRecord with quality scores
            output_file: Path to output filtered FASTA file

        Returns:
            Filtered SeqRecord if it meets length requirements, None otherwise
        """
        if "phred_quality" not in record.letter_annotations:
            raise ValueError("Record does not contain phred_quality annotations")

        qualities = record.letter_annotations["phred_quality"]
        sequence = str(record.seq)

        # Filter low quality bases
        filtered_seq = "".join(
            [
                base if qual >= self.min_quality else "N"
                for base, qual in zip(sequence, qualities)
            ]
        )

        # Check minimum sequence length (count only valid nucleotides)
        valid_bases = sum(1 for base in filtered_seq if base in "ATCG")
        if valid_bases < self.min_sequence_length:
            logger.warning(
                f"Sequence {record.id} too short after filtering: {valid_bases} valid bases (minimum: {self.min_sequence_length})"
            )
            return None

        # Create filtered record
        from Bio.Seq import Seq

        filtered_record = record[:]
        filtered_record.seq = Seq(filtered_seq)

        # Ensure output directory exists
        output_file.parent.mkdir(parents=True, exist_ok=True)

        # Write filtered FASTA
        SeqIO.write(filtered_record, output_file, "fasta")
        logger.info(
            f"Wrote filtered FASTA file: {output_file} ({valid_bases} valid bases)"
        )

        return filtered_record

    def generate_quality_plot(
        self,
        record: SeqRecord,
        output_file: Path,
        figure_size: Tuple[int, int] = (12, 4),
    ) -> None:
        """
        Generate quality score plot for the sequence.

        Args:
            record: SeqRecord with quality scores
            output_file: Path to output plot file
            figure_size: Figure size as (width, height)
        """
        if "phred_quality" not in record.letter_annotations:
            raise ValueError("Record does not contain phred_quality annotations")

        qualities = record.letter_annotations["phred_quality"]

        # Ensure output directory exists
        output_file.parent.mkdir(parents=True, exist_ok=True)

        # Create plot
        plt.figure(figsize=figure_size)
        plt.plot(qualities, marker=".", linestyle="-", markersize=1)
        plt.axhline(
            self.min_quality,
            color="red",
            linestyle="--",
            label=f"Quality threshold ({self.min_quality})",
        )
        plt.title(f"Quality scores for {record.id}")
        plt.xlabel("Base position")
        plt.ylabel("Phred Quality")
        plt.legend()
        plt.grid(True, alpha=0.3)
        plt.tight_layout()

        # Save plot
        plt.savefig(output_file, dpi=300, bbox_inches="tight")
        plt.close()
        logger.info(f"Generated quality plot: {output_file}")

    def process_ab1_file(
        self,
        ab1_file: Path,
        fasta_output: Path,
        filtered_output: Path,
        plot_output: Path,
    ) -> Tuple[SeqRecord, Optional[SeqRecord]]:
        """
        Complete processing of an AB1 file: convert, filter, and plot.

        Args:
            ab1_file: Path to input AB1 file
            fasta_output: Path to raw FASTA output
            filtered_output: Path to filtered FASTA output
            plot_output: Path to quality plot output

        Returns:
            Tuple of (raw_record, filtered_record or None if too short)
        """
        logger.info(f"Processing AB1 file: {ab1_file}")

        # Convert to FASTA
        raw_record = self.convert_to_fasta(ab1_file, fasta_output)

        # Filter by quality and check length
        filtered_record = self.filter_by_quality(raw_record, filtered_output)

        # Generate quality plot (always generate, even for short sequences)
        self.generate_quality_plot(raw_record, plot_output)

        if filtered_record is None:
            logger.warning(
                f"Sequence from {ab1_file} excluded due to insufficient length"
            )
        else:
            logger.info(f"Completed processing: {ab1_file}")

        return raw_record, filtered_record
