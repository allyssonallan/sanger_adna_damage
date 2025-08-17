#!/usr/bin/env python3
"""
Convert individual HVS consensus files to HSD format for HaploGrep analysis.
This version processes individual HVS consensus files rather than merged sequences.

Usage:
    python convert_hvs_consensus_to_hsd.py <consensus_directory> <output_hsd_file>
"""

import sys
import os
import re
from pathlib import Path
from Bio import SeqIO
from typing import Dict, List
from collections import defaultdict


class HVSConsensusToHSDConverter:
    """Convert individual HVS consensus files to HSD format."""
    
    def __init__(self, reference_file: str = "ref/rCRS.fasta"):
        """
        Initialize converter with reference sequence.
        
        Args:
            reference_file: Path to rCRS reference FASTA file
        """
        self.reference_file = reference_file
        self.reference_seq = ""
        self.hvs_regions = {
            'HVS1': {'start': 16024, 'end': 16365},
            'HVS2': {'start': 57, 'end': 372},
            'HVS3': {'start': 438, 'end': 574}
        }
        self.load_reference()
    
    def load_reference(self):
        """Load reference sequence from FASTA file."""
        try:
            with open(self.reference_file, 'r') as f:
                for record in SeqIO.parse(f, "fasta"):
                    self.reference_seq = str(record.seq)
                    break
            print(f"Reference sequence loaded: {len(self.reference_seq)} bases")
        except Exception as e:
            raise Exception(f"Error loading reference: {e}")
    
    def _clean_sequence(self, sequence: str) -> str:
        """Clean sequence by removing non-ATCG characters."""
        return re.sub(r'[^ATCG]', '', sequence.upper())
    
    def parse_sample_name(self, filename: str) -> tuple:
        """
        Parse sample name and HVS region from filename.
        
        Returns:
            (sample_id, hvs_region) tuple
        """
        # Remove extension and _consensus suffix
        basename = filename.replace('.fasta', '').replace('_consensus', '')
        
        # Extract HVS region
        hvs_region = None
        if 'HVS1' in basename.upper():
            hvs_region = 'HVS1'
        elif 'HVS2' in basename.upper():
            hvs_region = 'HVS2'
        elif 'HVS3' in basename.upper():
            hvs_region = 'HVS3'
        
        # Extract sample ID (everything before the HVS region)
        if hvs_region:
            pattern = f'_{hvs_region.lower()}'
            sample_id = basename.lower().replace(pattern, '')
        else:
            sample_id = basename.lower()
        
        return sample_id, hvs_region
    
    def find_variants(self, sample_seq: str, hvs_region: str) -> List[str]:
        """
        Find variants between sample and reference HVS region.
        
        Args:
            sample_seq: Cleaned sample sequence
            hvs_region: HVS region (HVS1, HVS2, or HVS3)
            
        Returns:
            List of variants in HSD format
        """
        variants = []
        
        if hvs_region not in self.hvs_regions:
            return variants
        
        # Get reference sequence for this HVS region
        region_info = self.hvs_regions[hvs_region]
        ref_start = region_info['start'] - 1  # Convert to 0-based
        ref_end = region_info['end']
        ref_region = self.reference_seq[ref_start:ref_end]
        ref_clean = self._clean_sequence(ref_region)
        
        # Compare sequences
        min_length = min(len(sample_seq), len(ref_clean))
        
        # Skip very short sequences
        if min_length < 50:
            return variants
        
        # Limit comparison to reasonable length to avoid excessive variants
        # Most ancient DNA HVS sequences should be similar to reference
        max_differences = min_length // 10  # Allow max 10% differences
        differences_found = 0
        
        for pos in range(min_length):
            if differences_found >= max_differences:
                break
                
            ref_base = ref_clean[pos]
            sample_base = sample_seq[pos]
            
            if ref_base != sample_base:
                # Convert back to mitochondrial coordinates
                mt_position = region_info['start'] + pos
                variant = f"{mt_position}{sample_base}"
                variants.append(variant)
                differences_found += 1
        
        return variants
    
    def process_consensus_directory(self, consensus_dir: str) -> Dict[str, List[str]]:
        """
        Process all consensus files in directory.
        
        Args:
            consensus_dir: Path to directory containing consensus FASTA files
            
        Returns:
            Dictionary mapping sample_id to list of variants
        """
        consensus_path = Path(consensus_dir)
        sample_variants = defaultdict(list)
        
        # Find all consensus files
        consensus_files = list(consensus_path.glob("*_consensus.fasta"))
        
        if not consensus_files:
            consensus_files = list(consensus_path.glob("*.fasta"))
        
        print(f"Found {len(consensus_files)} consensus files")
        
        for fasta_file in consensus_files:
            try:
                # Parse filename
                sample_id, hvs_region = self.parse_sample_name(fasta_file.name)
                
                if not hvs_region:
                    print(f"Skipping {fasta_file.name}: No HVS region detected")
                    continue
                
                # Read sequence
                with open(fasta_file, 'r') as f:
                    for record in SeqIO.parse(f, "fasta"):
                        sample_seq = self._clean_sequence(str(record.seq))
                        
                        if len(sample_seq) < 50:
                            print(f"Skipping {fasta_file.name}: Sequence too short ({len(sample_seq)}bp)")
                            continue
                        
                        # Find variants
                        variants = self.find_variants(sample_seq, hvs_region)
                        
                        if variants:
                            sample_variants[sample_id].extend(variants)
                            print(f"📁 {fasta_file.name}: {len(variants)} variants ({hvs_region})")
                        else:
                            print(f"📁 {fasta_file.name}: No variants found ({hvs_region})")
                        
                        break  # Only process first sequence in file
                        
            except Exception as e:
                print(f"Error processing {fasta_file.name}: {e}")
                continue
        
        return dict(sample_variants)
    
    def write_hsd_file(self, sample_variants: Dict[str, List[str]], output_file: str):
        """
        Write variants to HSD format file.
        
        Args:
            sample_variants: Dictionary mapping sample_id to variants
            output_file: Output HSD file path
        """
        with open(output_file, 'w') as f:
            for sample_id, variants in sample_variants.items():
                if variants:
                    # Sort variants by position
                    def extract_position(variant):
                        match = re.match(r'(\d+)', variant)
                        return int(match.group(1)) if match else 0
                    
                    sorted_variants = sorted(variants, key=extract_position)
                    
                    # Write HSD line
                    variants_str = '\t'.join(sorted_variants)
                    f.write(f"{sample_id}\t1-16569\t?\t{variants_str}\n")
                else:
                    # No variants found
                    f.write(f"{sample_id}\t1-16569\t?\n")
        
        print(f"✅ HSD file written: {output_file}")
        print(f"📋 Processed {len(sample_variants)} samples")


def main():
    if len(sys.argv) != 3:
        print("Usage: python convert_hvs_consensus_to_hsd.py <consensus_directory> <output_hsd_file>")
        sys.exit(1)
    
    consensus_dir = sys.argv[1]
    output_file = sys.argv[2]
    
    if not os.path.exists(consensus_dir):
        print(f"Error: Consensus directory not found: {consensus_dir}")
        sys.exit(1)
    
    print("🧬 Converting HVS consensus files to HSD format...")
    
    try:
        converter = HVSConsensusToHSDConverter()
        sample_variants = converter.process_consensus_directory(consensus_dir)
        converter.write_hsd_file(sample_variants, output_file)
        
        print("\n✅ Conversion complete!")
        print(f"📄 HSD file ready for HaploGrep analysis: {output_file}")
        print("🌐 HaploGrep web interface: https://haplogrep.i-med.ac.at/")
        
    except Exception as e:
        print(f"Error: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
