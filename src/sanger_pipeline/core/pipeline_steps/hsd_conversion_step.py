"""
HSD conversion step handler for the Sanger pipeline.

This module handles conversion of consensus sequences to HSD format for haplogroup analysis.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Dict

logger = logging.getLogger(__name__)


class HSDConversionStep:
    """Handles HSD format conversion for haplogroup analysis."""

    def __init__(self, config: Dict):
        """
        Initialize the HSD conversion step.

        Args:
            config: Pipeline configuration dictionary
        """
        self.config = config

    def execute(self, directories: Dict[str, Path]) -> Dict[str, int]:
        """
        Execute HSD conversion step.

        Args:
            directories: Pipeline output directories

        Returns:
            Dictionary with conversion statistics
        """
        logger.info("Step 6: Converting consensus sequences to HSD format")

        # Check if HSD conversion is enabled in config
        hsd_config = self.config.get("hsd_conversion", {})
        if not hsd_config.get("enabled", True):  # Default to enabled
            logger.info("HSD conversion disabled in configuration")
            return {"total_samples": 0, "samples_with_variants": 0, "total_variants": 0}

        consensus_dir = directories["consensus"]
        consensus_files = list(consensus_dir.glob("*_consensus.fasta"))

        if not consensus_files:
            logger.warning("No consensus files found for HSD conversion")
            return {"total_samples": 0, "samples_with_variants": 0, "total_variants": 0}

        try:
            # Initialize regional HSD converter
            reference_file = hsd_config.get("reference_file", "ref/rCRS.fasta")

            from ...utils.regional_hsd_converter import HybridRegionalHSDConverter

            converter = HybridRegionalHSDConverter(reference_file)

            # Process consensus files
            logger.info(
                f"Processing {len(consensus_files)} consensus files for HSD conversion"
            )
            sample_variants = converter.process_consensus_directory(str(consensus_dir))

            # Write HSD output
            hsd_output = directories["output"] / "haplogroup_analysis.hsd"
            converter.write_hsd_file(sample_variants, str(hsd_output))

            # Calculate statistics
            total_samples = len(sample_variants)
            samples_with_variants = sum(
                1 for variants in sample_variants.values() if variants
            )
            total_variants = sum(len(variants) for variants in sample_variants.values())

            # Log summary
            logger.info("HSD conversion completed:")
            logger.info(f"  - Processed {total_samples} samples")
            logger.info(f"  - {samples_with_variants} samples with variants")
            logger.info(f"  - {total_variants} total variants identified")
            logger.info(f"  - Output written to: {hsd_output}")

            return {
                "total_samples": total_samples,
                "samples_with_variants": samples_with_variants,
                "total_variants": total_variants,
            }

        except Exception as e:
            logger.error(f"HSD conversion failed: {e}")
            # Don't fail the entire pipeline for HSD conversion errors
            logger.warning("Continuing pipeline without HSD conversion")
            return {"total_samples": 0, "samples_with_variants": 0, "total_variants": 0}
