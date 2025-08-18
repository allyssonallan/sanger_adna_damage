"""
Consensus building step handler for the Sanger pipeline.

This module handles alignment and consensus generation for paired HVS region reads.
"""

import logging
from pathlib import Path
from typing import Dict

from ..consensus_builder import ConsensusBuilder

logger = logging.getLogger(__name__)


class ConsensusStep:
    """Handles alignment and consensus generation for paired reads."""

    def __init__(self, consensus_builder: ConsensusBuilder):
        """
        Initialize the consensus step.

        Args:
            consensus_builder: Configured ConsensusBuilder instance
        """
        self.consensus_builder = consensus_builder

    def execute(self, directories: Dict[str, Path]) -> Dict[str, int]:
        """
        Execute consensus generation step.

        Args:
            directories: Pipeline output directories

        Returns:
            Dictionary with consensus generation statistics
        """
        logger.info(
            "Step 2: Generating consensus sequences from paired reads by HVS region"
        )

        fasta_files = list(directories["fasta"].glob("*.fasta"))
        if not fasta_files:
            logger.warning("No FASTA files found for consensus generation")
            return {"processed_samples": 0, "total_consensus": 0}

        # Group files by sample and HVS region
        sample_hvs_groups = self._group_files_by_sample_and_region(fasta_files)

        # Process each sample's HVS regions
        processed_samples = 0
        total_consensus = 0

        for sample_name, hvs_regions in sample_hvs_groups.items():
            logger.info(
                f"Processing sample {sample_name} with regions: {list(hvs_regions.keys())}"
            )

            for hvs_region, files in hvs_regions.items():
                if "F" not in files or "R" not in files:
                    logger.warning(
                        f"Missing paired reads for {sample_name} {hvs_region} (F:{files.get('F')} R:{files.get('R')})"
                    )
                    continue

                forward_file = files["F"]
                reverse_file = files["R"]

                # Output files named with sample and HVS region
                alignment_output = (
                    directories["aligned"] / f"{sample_name}_{hvs_region}_aligned.fasta"
                )
                consensus_output = (
                    directories["consensus"]
                    / f"{sample_name}_{hvs_region}_consensus.fasta"
                )

                try:
                    self.consensus_builder.process_paired_reads(
                        forward_file,
                        reverse_file,
                        alignment_output,
                        consensus_output,
                        f"{sample_name}_{hvs_region}",
                    )
                    total_consensus += 1
                    logger.debug(f"Generated consensus for {sample_name} {hvs_region}")
                except Exception as e:
                    logger.error(
                        f"Failed to process paired reads for {sample_name} {hvs_region}: {e}"
                    )
                    continue

            if hvs_regions:  # If this sample had any HVS regions processed
                processed_samples += 1

        logger.info(
            f"Generated consensus for {total_consensus} HVS regions across {processed_samples} samples"
        )

        return {
            "processed_samples": processed_samples,
            "total_consensus": total_consensus,
        }

    def _group_files_by_sample_and_region(
        self, fasta_files: list
    ) -> Dict[str, Dict[str, Dict[str, Path]]]:
        """
        Group FASTA files by sample name and HVS region.

        Args:
            fasta_files: List of FASTA file paths

        Returns:
            Nested dictionary: {sample_name: {hvs_region: {direction: file_path}}}
        """
        sample_hvs_groups = {}

        for fasta_file in fasta_files:
            name = fasta_file.stem

            # Parse filename pattern: SAMPLE_NAME_HVS#-Direction
            if "_HVS" in name and name.endswith(("-F", "-R")):
                # Split at HVS region
                hvs_start = name.rfind("_HVS")
                sample_name = name[:hvs_start]
                hvs_part = name[hvs_start + 1 :]  # Remove the '_'

                # Extract HVS region and direction
                if hvs_part.endswith("-F"):
                    hvs_region = hvs_part[:-2]  # Remove '-F'
                    direction = "F"
                elif hvs_part.endswith("-R"):
                    hvs_region = hvs_part[:-2]  # Remove '-R'
                    direction = "R"
                else:
                    logger.warning(f"Unknown direction pattern in {name}")
                    continue

                # Initialize nested dictionary structure
                if sample_name not in sample_hvs_groups:
                    sample_hvs_groups[sample_name] = {}
                if hvs_region not in sample_hvs_groups[sample_name]:
                    sample_hvs_groups[sample_name][hvs_region] = {}

                sample_hvs_groups[sample_name][hvs_region][direction] = fasta_file
            else:
                logger.warning(f"Skipping file with unexpected naming pattern: {name}")

        return sample_hvs_groups
