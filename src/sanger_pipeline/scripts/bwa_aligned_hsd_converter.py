#!/usr/bin/env python3
"""
BWA-aligned HSD Converter for Sanger aDNA Pipeline.

This converter uses BWA for proper alignment of consensus sequences to the reference
before variant calling, which should provide much more reliable results than 
simple string matching approaches.
"""

import os
import sys
import subprocess
import tempfile
from pathlib import Path
from Bio import SeqIO
from typing import Dict, List
from collections import defaultdict
import re


class BWAAlignedHSDConverter:
    """Convert HVS consensus files to HSD format using BWA alignment."""
    
    def __init__(self, reference_file: str = "ref/rCRS.fasta"):
        """
        Initialize converter with reference sequence.
        
        Args:
            reference_file: Path to rCRS reference FASTA file
        """
        self.reference_file = reference_file
        self.reference_seq = ""
        
        # HVS regions with their coordinate ranges
        self.hvs_regions = {
            'HVS1': {'start': 16024, 'end': 16365},
            'HVS2': {'start': 57, 'end': 372},
            'HVS3': {'start': 438, 'end': 574}
        }
        
        # Primer sequences to remove from consensus sequences before alignment
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
        self.prepare_bwa_index()
    
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
    
    def prepare_bwa_index(self):
        """Prepare BWA index for reference sequence if it doesn't exist."""
        index_files = [
            f"{self.reference_file}.amb",
            f"{self.reference_file}.ann", 
            f"{self.reference_file}.bwt",
            f"{self.reference_file}.pac",
            f"{self.reference_file}.sa"
        ]
        
        # Check if index files exist
        if not all(os.path.exists(f) for f in index_files):
            print("🔨 Creating BWA index for reference...")
            try:
                subprocess.run([
                    'bwa', 'index', self.reference_file
                ], check=True, capture_output=True, text=True)
                print("✅ BWA index created successfully")
            except subprocess.CalledProcessError as e:
                raise Exception(f"Error creating BWA index: {e}")
        else:
            print("✅ BWA index found")
    
    def _remove_primers(self, sequence: str, hvs_region: str) -> str:
        """
        Remove primer sequences from consensus sequence before alignment.
        
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
    
    def align_sequence_with_bwa(self, sequence: str, sample_id: str) -> str:
        """
        Align sequence to reference using BWA.
        
        Args:
            sequence: Input sequence to align
            sample_id: Sample identifier for temporary files
            
        Returns:
            SAM alignment output
        """
        with tempfile.TemporaryDirectory() as temp_dir:
            # Write sequence to temporary FASTA file
            temp_fasta = os.path.join(temp_dir, f"{sample_id}.fasta")
            with open(temp_fasta, 'w') as f:
                f.write(f">{sample_id}\n{sequence}\n")
            
            # Run BWA alignment
            try:
                result = subprocess.run([
                    'bwa', 'mem', 
                    '-T', '20',  # Minimum score threshold
                    '-k', '15',  # Minimum seed length
                    '-r', '1.5', # Re-seeding trigger
                    self.reference_file,
                    temp_fasta
                ], capture_output=True, text=True, check=True)
                
                return result.stdout
                
            except subprocess.CalledProcessError as e:
                print(f"⚠️  BWA alignment failed for {sample_id}: {e}")
                return ""
    
    def parse_sam_for_variants(self, sam_output: str, hvs_region: str) -> List[str]:
        """
        Parse SAM output to extract variants in the specified HVS region.
        
        Args:
            sam_output: BWA SAM output
            hvs_region: HVS region (HVS1, HVS2, or HVS3)
            
        Returns:
            List of variants in HSD format
        """
        variants = []
        
        if hvs_region not in self.hvs_regions:
            return variants
        
        region_info = self.hvs_regions[hvs_region]
        region_start = region_info['start']
        region_end = region_info['end']
        
        for line in sam_output.split('\n'):
            if line.startswith('@') or not line.strip():
                continue
            
            parts = line.split('\t')
            if len(parts) < 11:
                continue
            
            # Parse SAM fields
            flag = int(parts[1])
            pos = int(parts[3])
            cigar = parts[5]
            seq = parts[9]
            
            # Skip unmapped reads
            if flag & 0x4:
                continue
            
            # Parse CIGAR string and extract variants
            variants.extend(self._extract_variants_from_cigar(
                pos, cigar, seq, region_start, region_end
            ))
        
        return sorted(set(variants))
    
    def _extract_variants_from_cigar(self, pos: int, cigar: str, seq: str, 
                                   region_start: int, region_end: int) -> List[str]:
        """
        Extract variants from CIGAR string within specified region.
        
        Args:
            pos: Alignment start position (1-based)
            cigar: CIGAR string
            seq: Query sequence
            region_start: Start of HVS region
            region_end: End of HVS region
            
        Returns:
            List of variants in the specified region
        """
        variants = []
        ref_pos = pos
        seq_pos = 0
        
        # Parse CIGAR operations
        cigar_ops = re.findall(r'(\d+)([MIDNSHPX=])', cigar)
        
        for length, operation in cigar_ops:
            length = int(length)
            
            if operation == 'M':  # Match/mismatch
                # Check for mismatches in this region
                for i in range(length):
                    current_ref_pos = ref_pos + i
                    current_seq_pos = seq_pos + i
                    
                    # Check if position is within our HVS region
                    if region_start <= current_ref_pos <= region_end:
                        if current_seq_pos < len(seq):
                            seq_base = seq[current_seq_pos]
                            ref_base = self.reference_seq[current_ref_pos - 1]  # Convert to 0-based
                            
                            # Only report clear differences
                            if (seq_base in 'ATCG' and ref_base in 'ATCG' and 
                                seq_base != ref_base):
                                variant = f"{current_ref_pos}{seq_base}"
                                variants.append(variant)
                
                ref_pos += length
                seq_pos += length
                
            elif operation == 'I':  # Insertion
                # Insertions relative to reference (use decimal notation)
                if region_start <= ref_pos <= region_end:
                    insertion_seq = seq[seq_pos:seq_pos + length]
                    if all(base in 'ATCG' for base in insertion_seq):
                        # Use decimal notation: position.1base, position.2base, etc.
                        for i, base in enumerate(insertion_seq):
                            variant = f"{ref_pos}.{i+1}{base}"
                            variants.append(variant)
                seq_pos += length
                
            elif operation == 'D':  # Deletion
                # Deletions from reference (skip for now as format unclear)
                ref_pos += length
                
            elif operation in 'SH':  # Soft/hard clipping
                if operation == 'S':
                    seq_pos += length
                # Hard clipping doesn't advance sequence position
                
            elif operation == 'N':  # Skipped region
                ref_pos += length
                
        return variants
    
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
    
    def process_consensus_directory(self, consensus_dir: str) -> Dict[str, List[str]]:
        """
        Process all consensus files in directory using BWA alignment.
        
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
                        raw_sequence = str(record.seq)
                        
                        # Remove primers BEFORE alignment (best practice)
                        sequence = self._remove_primers(raw_sequence, hvs_region)
                        
                        if len(sequence) < 50:
                            print(f"Skipping {fasta_file.name}: Sequence too short after primer removal ({len(sequence)}bp)")
                            continue
                        
                        # Align with BWA (primer-free sequence)
                        sam_output = self.align_sequence_with_bwa(sequence, f"{sample_id}_{hvs_region}")
                        
                        if not sam_output:
                            print(f"⚠️  No alignment for {fasta_file.name}")
                            continue
                        
                        # Extract variants from alignment
                        variants = self.parse_sam_for_variants(sam_output, hvs_region)
                        
                        if variants:
                            sample_variants[sample_id].extend(variants)
                            print(f"📁 {fasta_file.name}: {len(variants)} variants ({hvs_region}) - BWA aligned (primers removed)")
                        else:
                            print(f"📁 {fasta_file.name}: No variants found ({hvs_region}) - BWA aligned (primers removed)")
                        
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
            sample_regions: List of HVS regions this sample covers
            
        Returns:
            Range string in HSD format
        """
        ranges = []
        
        # Add ranges for each HVS region the sample has
        for region in sample_regions:
            if region in self.hvs_regions:
                region_info = self.hvs_regions[region]
                ranges.append(f"{region_info['start']}-{region_info['end']}")
        
        return ";".join(ranges) if ranges else "16024-16365"
    
    def write_hsd_file(self, sample_variants: Dict[str, List[str]], output_file: str):
        """
        Write variants to HSD format file.
        
        Args:
            sample_variants: Dictionary mapping sample_id to list of variants
            output_file: Path to output HSD file
        """
        with open(output_file, 'w') as f:
            for sample_id, variants in sample_variants.items():
                if not variants:
                    continue
                
                # Determine which regions this sample covers based on variant positions
                sample_regions = []
                for variant in variants:
                    # Extract position from variant (format: 16111T or 73G or 315.1C)
                    pos_match = re.match(r'(\d+)', variant)
                    if pos_match:
                        pos = int(pos_match.group(1))
                        for region, coords in self.hvs_regions.items():
                            if coords['start'] <= pos <= coords['end']:
                                if region not in sample_regions:
                                    sample_regions.append(region)
                
                # Get range string
                range_str = self.get_actual_range(variants, sample_regions)
                
                # Format variants (space-separated)
                variants_str = " ".join(sorted(variants))
                
                # Write in HSD format: SampleID\tRange\tHaplogroup\tVariants
                f.write(f"{sample_id}\t{range_str}\t?\t{variants_str}\n")
        
        print(f"✅ HSD file written: {output_file}")
        print(f"📋 Processed {len(sample_variants)} samples")


def main():
    """Main function to run BWA-aligned HSD conversion."""
    if len(sys.argv) != 3:
        print("Usage: python bwa_aligned_hsd_converter.py <consensus_dir> <output_hsd>")
        print("Example: python bwa_aligned_hsd_converter.py output_min30_q30/consensus/ output.hsd")
        sys.exit(1)
    
    consensus_dir = sys.argv[1]
    output_file = sys.argv[2]
    
    print("🧬 Converting HVS consensus files to HSD format using BWA alignment...")
    
    try:
        converter = BWAAlignedHSDConverter()
        sample_variants = converter.process_consensus_directory(consensus_dir)
        converter.write_hsd_file(sample_variants, output_file)
        
        print("\n✅ Conversion complete!")
        print(f"📄 HSD file ready for HaploGrep analysis: {output_file}")
        print("🌐 HaploGrep web interface: https://haplogrep.i-med.ac.at/")
        
    except Exception as e:
        print(f"❌ Error during conversion: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
