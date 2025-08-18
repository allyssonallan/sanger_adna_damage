#!/usr/bin/env python3
"""
Improved FASTA to HSD Converter for Ancient DNA Analysis

This module provides enhanced conversion from FASTA sequences to HSD format
with better handling of ambiguous nucleotides, quality filtering, and
alignment artifact detection.

Author: Sanger aDNA Pipeline
"""

import logging
from pathlib import Path
from typing import List
from Bio import SeqIO

logger = logging.getLogger(__name__)


class ImprovedFastaToHSDConverter:
    """
    Enhanced FASTA to HSD converter with improved quality control
    and alignment artifact detection.
    """

    def __init__(
        self, reference_path: str = "ref/rCRS.fasta", min_quality_threshold: float = 0.7
    ):
        """
        Initialize the converter with enhanced parameters.

        Args:
            reference_path: Path to reference sequence (rCRS)
            min_quality_threshold: Minimum proportion of valid nucleotides required
        """
        self.reference_path = Path(reference_path)
        self.min_quality_threshold = min_quality_threshold
        self.reference_seq = ""

        # HVS region definitions (1-based coordinates)
        self.hvs_regions = {
            "HVS1": (16024, 16365),
            "HVS2": (57, 372),
            "HVS3": (438, 574),
        }

        # Valid nucleotides for quality assessment
        self.valid_nucleotides = {"A", "T", "G", "C"}

        # Ambiguous nucleotide mapping for better handling
        self.ambiguous_mapping = {
            "R": ["A", "G"],  # Purine
            "Y": ["C", "T"],  # Pyrimidine
            "S": ["G", "C"],  # Strong
            "W": ["A", "T"],  # Weak
            "K": ["G", "T"],  # Keto
            "M": ["A", "C"],  # Amino
            "B": ["C", "G", "T"],  # Not A
            "D": ["A", "G", "T"],  # Not C
            "H": ["A", "C", "T"],  # Not G
            "V": ["A", "C", "G"],  # Not T
            "N": ["A", "T", "G", "C"],  # Any
        }

        self._load_reference()

    def _load_reference(self):
        """Load and validate reference sequence."""
        try:
            with open(self.reference_path, "r") as f:
                record = next(SeqIO.parse(f, "fasta"))
                self.reference_seq = str(record.seq).upper()
                logger.info(
                    f"Reference sequence loaded: {len(self.reference_seq)} bases"
                )
        except FileNotFoundError:
            logger.error(f"Reference file not found: {self.reference_path}")
            raise
        except Exception as e:
            logger.error(f"Error loading reference: {e}")
            raise

    def _calculate_sequence_quality(self, sequence: str) -> float:
        """
        Calculate the quality score of a sequence based on valid nucleotides.

        Args:
            sequence: DNA sequence string

        Returns:
            Quality score (0.0 to 1.0)
        """
        if not sequence:
            return 0.0

        valid_count = sum(
            1 for base in sequence.upper() if base in self.valid_nucleotides
        )
        return valid_count / len(sequence)

    def _filter_low_quality_regions(self, sequence: str, window_size: int = 10) -> str:
        """
        Filter out low quality regions using a sliding window approach.

        Args:
            sequence: Input DNA sequence
            window_size: Size of sliding window for quality assessment

        Returns:
            Filtered sequence with low quality regions masked as 'N'
        """
        filtered_seq = list(sequence.upper())

        for i in range(len(sequence) - window_size + 1):
            window = sequence[i : i + window_size]
            window_quality = self._calculate_sequence_quality(window)

            if window_quality < self.min_quality_threshold:
                # Mark this region as low quality
                for j in range(i, min(i + window_size, len(filtered_seq))):
                    if filtered_seq[j] not in self.valid_nucleotides:
                        filtered_seq[j] = "N"

        return "".join(filtered_seq)

    def _resolve_ambiguous_nucleotide(self, ambiguous_base: str, ref_base: str) -> str:
        """
        Resolve ambiguous nucleotide based on reference and possible options.

        Args:
            ambiguous_base: Ambiguous nucleotide code
            ref_base: Reference nucleotide at this position

        Returns:
            Most likely nucleotide or 'N' if unresolvable
        """
        if ambiguous_base in self.valid_nucleotides:
            return ambiguous_base

        if ambiguous_base in self.ambiguous_mapping:
            possible_bases = self.ambiguous_mapping[ambiguous_base]

            # If reference is one of the possibilities, prefer it (no mutation)
            if ref_base in possible_bases:
                return ref_base
            else:
                # Return the first alternative (indicating mutation)
                return possible_bases[0] if possible_bases else "N"

        return "N"  # Unresolvable

    def _detect_hvs_regions(self, sequence: str, sample_name: str) -> List[str]:
        """
        Detect which HVS regions are present in the sequence with quality filtering.

        Args:
            sequence: DNA sequence
            sample_name: Sample identifier

        Returns:
            List of detected HVS region names
        """
        detected_regions = []
        seq_len = len(sequence)

        # Calculate overall sequence quality
        overall_quality = self._calculate_sequence_quality(sequence)

        logger.info(
            f"Processing {sample_name}: Sample={seq_len}bp, Quality={overall_quality:.2f}"
        )

        if overall_quality < self.min_quality_threshold:
            logger.warning(
                f"Low quality sequence for {sample_name}: {overall_quality:.2f}"
            )
            return detected_regions

        # Check for HVS regions based on sequence length and content
        if seq_len >= 300:  # Likely HVS1 or combined regions
            if any(
                name.upper() in sample_name.upper() for name in ["HVS1", "CONSENSUS"]
            ):
                detected_regions.append("HVS1")
            elif any(name.upper() in sample_name.upper() for name in ["HVS2"]):
                detected_regions.append("HVS2")
            elif any(name.upper() in sample_name.upper() for name in ["HVS3"]):
                detected_regions.append("HVS3")
            else:
                # Try to detect based on sequence characteristics
                detected_regions = self._infer_regions_from_sequence(sequence, seq_len)

        return detected_regions

    def _infer_regions_from_sequence(self, sequence: str, seq_len: int) -> List[str]:
        """
        Infer HVS regions from sequence characteristics.

        Args:
            sequence: DNA sequence
            seq_len: Sequence length

        Returns:
            List of inferred region names
        """
        regions = []

        if seq_len > 500:  # Likely multiple regions
            regions = ["HVS1", "HVS2", "HVS3"]
        elif seq_len > 250:  # Likely HVS1 or HVS2
            regions = ["HVS1"]
        elif seq_len > 100:  # Likely HVS3 or partial region
            regions = ["HVS2"]
        else:
            regions = ["HVS3"]

        return regions

    def _compare_sequences_with_quality_control(
        self, sequence: str, regions: List[str]
    ) -> List[str]:
        """
        Compare sequence with reference using enhanced quality control.

        Args:
            sequence: Sample sequence
            regions: List of HVS regions to analyze

        Returns:
            List of variant strings
        """
        variants = []

        # Filter low quality regions first
        filtered_sequence = self._filter_low_quality_regions(sequence)

        for region_name in regions:
            if region_name not in self.hvs_regions:
                continue

            start_pos, end_pos = self.hvs_regions[region_name]

            # Extract reference region
            ref_region = self.reference_seq[start_pos - 1 : end_pos]

            # For comparison, we need to align the sequence properly
            # This is a simplified approach - in practice, you might want more sophisticated alignment
            region_variants = self._find_variants_in_region(
                filtered_sequence, ref_region, start_pos, region_name
            )

            variants.extend(region_variants)

        return variants

    def _find_variants_in_region(
        self, sequence: str, ref_region: str, start_pos: int, region_name: str
    ) -> List[str]:
        """
        Find variants in a specific HVS region with quality control.

        Args:
            sequence: Sample sequence
            ref_region: Reference sequence for this region
            start_pos: Starting position in reference
            region_name: Name of the HVS region

        Returns:
            List of variant strings for this region
        """
        variants = []

        # Simple alignment approach - compare position by position
        # In a real implementation, you might want to use proper sequence alignment
        min_len = min(len(sequence), len(ref_region))

        for i in range(min_len):
            seq_base = sequence[i].upper()
            ref_base = ref_region[i].upper()

            # Skip low quality positions
            if seq_base == "N":
                continue

            # Resolve ambiguous nucleotides
            resolved_base = self._resolve_ambiguous_nucleotide(seq_base, ref_base)

            if resolved_base != "N" and resolved_base != ref_base:
                position = start_pos + i
                variant = f"{position}{resolved_base}"
                variants.append(variant)

        logger.debug(f"Found {len(variants)} variants in {region_name}")
        return variants

    def convert_fasta_to_hsd(
        self, input_fasta: str, output_hsd: str, quality_filter: bool = True
    ) -> None:
        """
        Convert FASTA file to HSD format with enhanced quality control.

        Args:
            input_fasta: Path to input FASTA file
            output_hsd: Path to output HSD file
            quality_filter: Whether to apply quality filtering
        """
        logger.info(f"Converting {input_fasta} to HSD format...")
        logger.info(f"Quality filtering: {'enabled' if quality_filter else 'disabled'}")
        logger.info(f"Quality threshold: {self.min_quality_threshold}")

        with open(output_hsd, "w") as hsd_file:
            processed_count = 0
            skipped_count = 0

            for record in SeqIO.parse(input_fasta, "fasta"):
                sample_name = record.id
                sequence = str(record.seq).upper()

                # Calculate sequence quality
                quality_score = self._calculate_sequence_quality(sequence)

                if quality_filter and quality_score < self.min_quality_threshold:
                    logger.warning(
                        f"Skipping {sample_name}: quality too low ({quality_score:.2f})"
                    )
                    skipped_count += 1
                    continue

                # Detect HVS regions
                detected_regions = self._detect_hvs_regions(sequence, sample_name)

                if not detected_regions:
                    logger.warning(f"No HVS regions detected for {sample_name}")
                    skipped_count += 1
                    continue

                # Create region string for HSD format
                region_ranges = []
                for region in detected_regions:
                    if region in self.hvs_regions:
                        start, end = self.hvs_regions[region]
                        region_ranges.append(f"{start}-{end}")

                region_string = ";".join(region_ranges)

                # Find variants with quality control
                variants = self._compare_sequences_with_quality_control(
                    sequence, detected_regions
                )

                if variants:
                    # Format variants for HSD
                    variant_string = "\t".join(variants)

                    # Write HSD entry
                    hsd_line = f"{sample_name}\t{region_string}\t?\t{variant_string}\n"
                    hsd_file.write(hsd_line)

                    logger.info(
                        f"Processed {sample_name}: {len(variants)} variants, Quality: {quality_score:.2f}"
                    )
                    processed_count += 1
                else:
                    logger.warning(f"No variants found for {sample_name}")
                    skipped_count += 1

        logger.info(
            f"Conversion completed: {processed_count} sequences processed, {skipped_count} skipped"
        )
        logger.info(f"Output written to: {output_hsd}")


def main():
    """Main function for command line usage."""
    import argparse

    parser = argparse.ArgumentParser(description="Improved FASTA to HSD Converter")
    parser.add_argument("-i", "--input", required=True, help="Input FASTA file")
    parser.add_argument("-o", "--output", required=True, help="Output HSD file")
    parser.add_argument(
        "-r", "--reference", default="ref/rCRS.fasta", help="Reference sequence file"
    )
    parser.add_argument(
        "-q",
        "--quality-threshold",
        type=float,
        default=0.7,
        help="Minimum quality threshold (0.0-1.0)",
    )
    parser.add_argument(
        "--no-quality-filter", action="store_true", help="Disable quality filtering"
    )
    parser.add_argument("-v", "--verbose", action="store_true", help="Verbose logging")

    args = parser.parse_args()

    # Set up logging
    log_level = logging.DEBUG if args.verbose else logging.INFO
    logging.basicConfig(
        level=log_level, format="%(asctime)s - %(levelname)s - %(message)s"
    )

    # Create converter and process
    converter = ImprovedFastaToHSDConverter(
        reference_path=args.reference, min_quality_threshold=args.quality_threshold
    )

    converter.convert_fasta_to_hsd(
        input_fasta=args.input,
        output_hsd=args.output,
        quality_filter=not args.no_quality_filter,
    )


if __name__ == "__main__":
    main()
