"""
Region merging step handler for the Sanger pipeline.

This module handles merging of available HVS regions into final consensus sequences.
"""

import logging
from pathlib import Path
from typing import Dict

logger = logging.getLogger(__name__)


class RegionMergingStep:
    """Handles merging of HVS regions into final consensus sequences."""
    
    def execute(self, directories: Dict[str, Path]) -> Dict[str, int]:
        """
        Execute region merging step.
        
        Args:
            directories: Pipeline output directories
            
        Returns:
            Dictionary with merging statistics
        """
        logger.info("Step 3: Merging available HVS regions into final consensus sequences")

        consensus_files = list(directories["consensus"].glob("*_consensus.fasta"))
        if not consensus_files:
            logger.warning("No consensus files found for merging")
            return {"processed_samples": 0}

        # Group consensus files by sample
        sample_consensus = self._group_consensus_by_sample(consensus_files)

        # Process each sample
        processed_samples = 0
        for sample_name, hvs_files in sample_consensus.items():
            available_regions = sorted(hvs_files.keys())
            logger.info(f"Sample {sample_name} has regions: {available_regions}")
            
            try:
                # Read all available HVS consensus sequences
                sequences = self._read_consensus_sequences(hvs_files)
                
                # Create merged sequence name and content based on available regions
                merged_name, merged_sequence = self._create_merged_sequence(
                    sample_name, available_regions, sequences
                )
                
                # Write merged consensus
                self._write_merged_consensus(
                    directories["final"], merged_name, merged_sequence
                )
                
                processed_samples += 1
                logger.debug(f"Created merged sequence for {sample_name}: {merged_name}")
                
            except Exception as e:
                logger.error(f"Failed to merge regions for {sample_name}: {e}")
                continue

        logger.info(f"Created merged sequences for {processed_samples} samples")
        
        return {"processed_samples": processed_samples}
    
    def _group_consensus_by_sample(self, consensus_files: list) -> Dict[str, Dict[str, Path]]:
        """
        Group consensus files by sample name.
        
        Args:
            consensus_files: List of consensus file paths
            
        Returns:
            Dictionary mapping sample names to their HVS region files
        """
        sample_consensus = {}
        
        for consensus_file in consensus_files:
            name = consensus_file.stem.replace("_consensus", "")
            
            # Parse sample name and HVS region
            if "_HVS" in name:
                hvs_start = name.rfind("_HVS")
                sample_name = name[:hvs_start]
                hvs_region = name[hvs_start+1:]  # Remove the '_'
                
                if sample_name not in sample_consensus:
                    sample_consensus[sample_name] = {}
                
                sample_consensus[sample_name][hvs_region] = consensus_file
            else:
                logger.warning(f"Unexpected consensus file naming: {name}")
        
        return sample_consensus
    
    def _read_consensus_sequences(self, hvs_files: Dict[str, Path]) -> Dict[str, str]:
        """
        Read consensus sequences from files.
        
        Args:
            hvs_files: Dictionary mapping HVS regions to file paths
            
        Returns:
            Dictionary mapping HVS regions to sequences
        """
        sequences = {}
        
        for hvs_region, consensus_file in hvs_files.items():
            with open(consensus_file, 'r') as f:
                # Read the sequence content (skip FASTA header)
                lines = f.readlines()
                sequence_lines = [line.strip() for line in lines if not line.startswith('>')]
                sequences[hvs_region] = ''.join(sequence_lines)
        
        return sequences
    
    def _create_merged_sequence(self, sample_name: str, available_regions: list, 
                              sequences: Dict[str, str]) -> tuple:
        """
        Create merged sequence name and content.
        
        Args:
            sample_name: Name of the sample
            available_regions: List of available HVS regions
            sequences: Dictionary mapping regions to sequences
            
        Returns:
            Tuple of (merged_name, merged_sequence)
        """
        if len(available_regions) == 1:
            # Single region
            region_name = available_regions[0]
            merged_name = f"{sample_name}_{region_name}"
            merged_sequence = sequences[available_regions[0]]
        else:
            # Multiple regions - concatenate in order (HVS1, HVS2, HVS3)
            ordered_regions = []
            merged_sequences = []
            
            for hvs_num in ['HVS1', 'HVS2', 'HVS3']:
                if hvs_num in available_regions:
                    ordered_regions.append(hvs_num)
                    merged_sequences.append(sequences[hvs_num])
            
            region_name = "_".join(ordered_regions)
            merged_name = f"{sample_name}_{region_name}"
            merged_sequence = ''.join(merged_sequences)  # Concatenate sequences
        
        return merged_name, merged_sequence
    
    def _write_merged_consensus(self, final_dir: Path, merged_name: str, 
                              merged_sequence: str) -> None:
        """
        Write merged consensus sequence to file.
        
        Args:
            final_dir: Directory for final output files
            merged_name: Name for merged sequence
            merged_sequence: Merged sequence content
        """
        merged_output = final_dir / f"{merged_name}_merged.fasta"
        
        with open(merged_output, 'w') as f:
            f.write(f">{merged_name}_consensus\n")
            # Write sequence with line breaks every 80 characters
            for i in range(0, len(merged_sequence), 80):
                f.write(merged_sequence[i:i+80] + '\n')
