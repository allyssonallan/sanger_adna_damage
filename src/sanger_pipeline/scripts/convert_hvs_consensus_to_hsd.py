#!/usr/bin/env python3
"""
Convert individual HVS consensus files to HSD format for HaploGrep analysis.
This version processes individual HVS consensus files rather than merged sequences.

Usage:
    python -m sanger_pipeline.scripts.convert_hvs_consensus_to_hsd <consensus_directory> <output_hsd_file>
"""

import sys
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
        
        # Primer sequences to remove from consensus sequences
        self.primers = {
            'HVS1': {
                'forward': 'CACCATTAGCACCCAAAGCT',
                'reverse': 'CACCATCCTCCGTGAAATCA'  # Reverse complement of HVS1-R
            },
            'HVS2': {
                'forward': 'GGTCTATCACCCTATTAACCAC', 
                'reverse': 'TGGCGGTATGCACTTTTAACAG'  # Reverse complement of HVS2-R
            },
            'HVS3': {
                'forward': 'CCGCTTCTGGCCACAGCACT',
                'reverse': 'GTTTAGACGGGCTCACATCACC'  # Reverse complement of HVS3-R
            }
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
        """Clean sequence by converting to uppercase but keeping all characters for alignment."""
        return sequence.upper()
    
    def _remove_primers(self, sequence: str, hvs_region: str) -> str:
        """
        Remove primer sequences from consensus sequence.
        
        Args:
            sequence: Raw consensus sequence
            hvs_region: HVS region (HVS1, HVS2, or HVS3)
            
        Returns:
            Sequence with primers removed
        """
        if hvs_region not in self.primers:
            return sequence
        
        seq_upper = sequence.upper()
        forward_primer = self.primers[hvs_region]['forward']
        reverse_primer = self.primers[hvs_region]['reverse']
        
        # Remove forward primer from the beginning
        if seq_upper.startswith(forward_primer):
            seq_upper = seq_upper[len(forward_primer):]
        else:
            # Try to find primer with some mismatches (up to 2)
            for i in range(min(len(forward_primer) + 5, len(seq_upper))):
                window = seq_upper[i:i + len(forward_primer)]
                if len(window) == len(forward_primer):
                    matches = sum(1 for a, b in zip(window, forward_primer) if a == b)
                    if matches >= len(forward_primer) - 2:  # Allow up to 2 mismatches
                        seq_upper = seq_upper[i + len(forward_primer):]
                        break
        
        # Remove reverse primer from the end
        if seq_upper.endswith(reverse_primer):
            seq_upper = seq_upper[:-len(reverse_primer)]
        else:
            # Try to find primer with some mismatches (up to 2)
            for i in range(min(len(reverse_primer) + 5, len(seq_upper))):
                start_pos = len(seq_upper) - len(reverse_primer) - i
                if start_pos >= 0:
                    window = seq_upper[start_pos:start_pos + len(reverse_primer)]
                    if len(window) == len(reverse_primer):
                        matches = sum(1 for a, b in zip(window, reverse_primer) if a == b)
                        if matches >= len(reverse_primer) - 2:  # Allow up to 2 mismatches
                            seq_upper = seq_upper[:start_pos]
                            break
        
        return seq_upper
    
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
        Find variants between sample and reference HVS region using alignment.
        
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
        ref_clean = ref_region.upper()
        
        # Find the best alignment position in the sample sequence
        best_alignment_pos = self._find_alignment_position(sample_seq, ref_clean)
        
        if best_alignment_pos == -1:
            print(f"Warning: Could not align {hvs_region} region in sample")
            return variants
        
        # Extract the aligned portion of the sample sequence
        aligned_sample = sample_seq[best_alignment_pos:best_alignment_pos + len(ref_clean)]
        
        # Compare aligned sequences
        min_length = min(len(aligned_sample), len(ref_clean))
        
        for pos in range(min_length):
            ref_base = ref_clean[pos]
            sample_base = aligned_sample[pos]
            
            # Only report actual differences (variants)
            # Skip positions with ambiguous bases in sample
            if (sample_base in 'ATCG' and 
                ref_base in 'ATCG' and 
                ref_base != sample_base):
                
                # Convert back to mitochondrial coordinates
                mt_position = region_info['start'] + pos
                variant = f"{mt_position}{sample_base}"
                variants.append(variant)
        
        return variants
    
    def _find_alignment_position(self, sample_seq: str, ref_seq: str) -> int:
        """
        Find the best alignment position of reference sequence within sample sequence.
        
        Args:
            sample_seq: Sample sequence to search in
            ref_seq: Reference sequence to find
            
        Returns:
            Best alignment position, or -1 if no good alignment found
        """
        best_pos = -1
        best_score = 0
        min_score = 20  # Minimum number of matching bases required
        
        # Use a sliding window to find the best match
        window_size = min(50, len(ref_seq) // 3)  # Use first third of reference for alignment
        ref_window = ref_seq[:window_size]
        
        for start_pos in range(len(sample_seq) - window_size + 1):
            sample_window = sample_seq[start_pos:start_pos + window_size]
            
            # Count matches (only clear bases)
            matches = 0
            total_clear = 0
            
            for i in range(window_size):
                if (i < len(ref_window) and i < len(sample_window) and
                    sample_window[i] in 'ATCG' and ref_window[i] in 'ATCG'):
                    total_clear += 1
                    if sample_window[i] == ref_window[i]:
                        matches += 1
            
            # Calculate score as percentage of matches among clear bases
            if total_clear > 0:
                score = matches
                if score > best_score and score >= min_score:
                    best_score = score
                    best_pos = start_pos
        
        return best_pos
    
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
                        raw_seq = str(record.seq)
                        
                        # Remove primers first
                        seq_no_primers = self._remove_primers(raw_seq, hvs_region)
                        sample_seq = seq_no_primers.upper()
                        
                        if len(sample_seq) < 50:
                            print(f"Skipping {fasta_file.name}: Sequence too short after primer removal ({len(sample_seq)}bp)")
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
    
    def get_actual_range(self, sample_variants: List[str], sample_regions: List[str]) -> str:
        """
        Get actual range covered by sample based on regions and variants.
        
        Args:
            sample_variants: List of variants for this sample
            sample_regions: List of HVS regions for this sample
            
        Returns:
            Range string for HSD format
        """
        ranges = []
        
        # Add ranges for each HVS region the sample has
        for region in sample_regions:
            if region in self.hvs_regions:
                region_info = self.hvs_regions[region]
                ranges.append(f"{region_info['start']}-{region_info['end']}")
        
        # If no regions detected, use variants to determine coverage
        if not ranges and sample_variants:
            positions = []
            for variant in sample_variants:
                match = re.match(r'(\d+)', variant)
                if match:
                    positions.append(int(match.group(1)))
            
            if positions:
                min_pos = min(positions)
                max_pos = max(positions)
                # Determine which HVS regions are covered
                for region, coords in self.hvs_regions.items():
                    if min_pos >= coords['start'] and max_pos <= coords['end']:
                        ranges.append(f"{coords['start']}-{coords['end']}")
                        break
        
        return ";".join(ranges) if ranges else "16024-16365"  # Default to HVS1 if nothing detected

    def write_hsd_file(self, sample_variants: Dict[str, List[str]], output_file: str):
        """
        Write variants to HSD format file.
        
        Args:
            sample_variants: Dictionary mapping sample_id to variants
            output_file: Output HSD file path
        """
        # First, we need to collect region info for each sample
        sample_regions = defaultdict(list)
        
        # Re-scan files to get region information
        consensus_path = Path(output_file).parent.parent / "consensus"
        if consensus_path.exists():
            for fasta_file in consensus_path.glob("*_consensus.fasta"):
                sample_id, hvs_region = self.parse_sample_name(fasta_file.name)
                if hvs_region:
                    sample_regions[sample_id].append(hvs_region)
        
        with open(output_file, 'w') as f:
            for sample_id, variants in sample_variants.items():
                # Get actual range for this sample
                actual_range = self.get_actual_range(variants, sample_regions[sample_id])
                
                if variants:
                    # Sort variants by position
                    def extract_position(variant):
                        match = re.match(r'(\d+)', variant)
                        return int(match.group(1)) if match else 0
                    
                    sorted_variants = sorted(variants, key=extract_position)
                    
                    # Write HSD line with actual regional ranges
                    variants_str = ' '.join(sorted_variants)  # Space-separated variants
                    f.write(f"{sample_id}\t{actual_range}\t?\t{variants_str}\n")
                else:
                    # No variants found - matches reference perfectly
                    f.write(f"{sample_id}\t{actual_range}\t?\n")
        
        print(f"✅ HSD file written: {output_file}")
        print(f"📋 Processed {len(sample_variants)} samples")


def main():
    if len(sys.argv) != 3:
        print("Usage: python convert_hvs_consensus_to_hsd.py <consensus_directory> <output_hsd_file>")
        sys.exit(1)
    
    consensus_dir = sys.argv[1]
    output_file = sys.argv[2]
    
    if not Path(consensus_dir).exists():
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
