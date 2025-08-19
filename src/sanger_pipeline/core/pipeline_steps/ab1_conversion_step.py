"""
AB1 conversion step handler for the Sanger pipeline.

This module handles the conversion of AB1 files to FASTA format with quality filtering.
"""

import logging
from pathlib import Path
from typing import Dict

from ..enhanced_ab1_converter_fixed import EnhancedAB1Converter as AB1Converter

logger = logging.getLogger(__name__)


class AB1ConversionStep:
    """Handles AB1 file conversion with quality filtering."""

    def __init__(self, ab1_converter: AB1Converter):
        """
        Initialize the AB1 conversion step.

        Args:
            ab1_converter: Configured AB1Converter instance
        """
        self.ab1_converter = ab1_converter

    def execute(self, input_dir: Path, directories: Dict[str, Path]) -> Dict[str, int]:
        """
        Execute AB1 conversion step.

        Args:
            input_dir: Directory containing AB1 files
            directories: Pipeline output directories

        Returns:
            Dictionary with conversion statistics
        """
        logger.info(
            "Step 1: Converting AB1 files to FASTA with quality filtering and length validation"
        )

        ab1_files = list(input_dir.glob("*.ab1"))
        if not ab1_files:
            logger.warning(f"No AB1 files found in {input_dir}")
            return {"processed_files": 0, "filtered_files": 0, "excluded_files": 0}

        processed_files = 0
        filtered_files = 0
        excluded_files = 0

        for ab1_file in ab1_files:
            fasta_output = directories["fasta"] / f"{ab1_file.stem}.fasta"
            filtered_output = (
                directories["filtered"] / f"{ab1_file.stem}_filtered.fasta"
            )
            plot_output = directories["plots"] / f"{ab1_file.stem}_quality.png"

            try:
                # Process AB1 file: convert, filter, and plot using enhanced method
                raw_record, filtered_record, stats = self.ab1_converter.process_ab1_file_enhanced(
                    ab1_file, fasta_output, filtered_output, plot_output
                )
                processed_files += 1

                if filtered_record is not None:
                    filtered_files += 1
                else:
                    excluded_files += 1
                    # Remove the empty filtered file if it was created but sequence was too short
                    if filtered_output.exists():
                        filtered_output.unlink()

            except Exception as e:
                logger.error(f"Failed to process {ab1_file}: {e}")
                continue

        logger.info(
            f"Processed {processed_files} AB1 files: {filtered_files} passed filtering, {excluded_files} excluded (too short)"
        )

        return {
            "processed_files": processed_files,
            "filtered_files": filtered_files,
            "excluded_files": excluded_files,
        }
