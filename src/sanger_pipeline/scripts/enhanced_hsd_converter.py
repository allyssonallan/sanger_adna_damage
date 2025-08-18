#!/usr/bin/env python3
"""
Enhanced HSD converter with reference-aware alignment for improved variant calling.

This module provides both the original direct comparison method and an enhanced
reference-aligned method for converting consensus sequences to HSD format.

Usage:
    python -m sanger_pipeline.scripts.enhanced_hsd_converter <consensus_directory> <output_hsd_file> [--method=aligned]
"""

import sys
import re
import tempfile
import subprocess
from pathlib import Path
from Bio import SeqIO, AlignIO
from typing import Dict, List, Tuple, Optional
from collections import defaultdict


class EnhancedHSDConverter:
    """
    Enhanced HSD converter with reference alignment capabilities.

    Provides two methods:
    1. Direct comparison (original method)
    2. Reference-aligned comparison (enhanced method)
    """

    def __init__(
        self, reference_file: str = "ref/rCRS.fasta", alignment_tool: str = "mafft"
    ):
        """
        Initialize enhanced converter.

        Args:
            reference_file: Path to rCRS reference FASTA file
            alignment_tool: Alignment tool to use (default: mafft)
        """
        self.reference_file = reference_file
        self.alignment_tool = alignment_tool
        self.reference_seq = ""
        self.hvs_regions = {
            "HVS1": {"start": 16024, "end": 16365},
            "HVS2": {"start": 57, "end": 372},
            "HVS3": {"start": 438, "end": 574},
        }
        self.load_reference()

    def load_reference(self):
        """Load the reference sequence."""
        try:
            reference_path = Path(self.reference_file)
            if not reference_path.exists():
                # Try relative to project root
                reference_path = (
                    Path(__file__).parent.parent.parent.parent / self.reference_file
                )

            if reference_path.exists():
                record = SeqIO.read(reference_path, "fasta")
                self.reference_seq = str(record.seq).upper()
                print(f"Loaded reference sequence: {len(self.reference_seq)} bp")
            else:
                raise FileNotFoundError(
                    f"Reference file not found: {self.reference_file}"
                )
        except Exception as e:
            print(f"Error loading reference: {e}")
            sys.exit(1)

    def load_consensus_sequences(self, consensus_dir: str) -> Dict[str, Dict[str, str]]:
        """
        Load consensus sequences from directory.

        Args:
            consensus_dir: Directory containing consensus FASTA files

        Returns:
            Dictionary mapping sample_id -> region -> sequence
        """
        sequences = defaultdict(dict)
        consensus_path = Path(consensus_dir)

        if not consensus_path.exists():
            raise FileNotFoundError(f"Consensus directory not found: {consensus_dir}")

        fasta_files = list(consensus_path.glob("*.fasta")) + list(
            consensus_path.glob("*.fa")
        )
        print(f"Found {len(fasta_files)} FASTA files in {consensus_dir}")

        for fasta_file in fasta_files:
            try:
                # Parse filename to extract sample ID and region
                filename = fasta_file.stem
                sample_id, region = self.parse_filename(filename)

                if sample_id and region:
                    record = SeqIO.read(fasta_file, "fasta")
                    sequences[sample_id][region] = str(record.seq).upper()
                    print(f"Loaded {sample_id} {region}: {len(record.seq)} bp")
                else:
                    print(f"Warning: Could not parse filename {filename}")

            except Exception as e:
                print(f"Error reading {fasta_file}: {e}")
                continue

        return dict(sequences)

    def parse_filename(self, filename: str) -> Tuple[Optional[str], Optional[str]]:
        """
        Parse filename to extract sample ID and HVS region.

        Args:
            filename: Filename without extension

        Returns:
            Tuple of (sample_id, region) or (None, None) if parsing fails
        """
        # Common patterns for HVS consensus files
        patterns = [
            r"(.+)_HVS([123])_consensus",  # sample_HVS1_consensus
            r"(.+)_hvs([123])_consensus",  # sample_hvs1_consensus
            r"(.+)_HVS([123])",  # sample_HVS1
            r"(.+)_hvs([123])",  # sample_hvs1
            r"(.+)-(HVS[123])",  # sample-HVS1
            r"(.+)-(hvs[123])",  # sample-hvs1
        ]

        for pattern in patterns:
            match = re.match(pattern, filename)
            if match:
                sample_id = match.group(1)
                region_num = match.group(2)
                if region_num.upper() in ["1", "2", "3"]:
                    region = f"HVS{region_num.upper()}"
                elif region_num.upper() in ["HVS1", "HVS2", "HVS3"]:
                    region = region_num.upper()
                else:
                    continue
                return sample_id, region

        return None, None

    def create_alignment(self, sequences: List[str], temp_dir: str) -> Optional[str]:
        """
        Create multiple sequence alignment using external tool.

        Args:
            sequences: List of sequences to align
            temp_dir: Temporary directory for files

        Returns:
            Path to alignment file or None if failed
        """
        # Create temporary input file
        input_file = Path(temp_dir) / "input.fasta"
        with open(input_file, "w") as f:
            for i, seq in enumerate(sequences):
                f.write(f">seq_{i}\n{seq}\n")

        # Create alignment
        output_file = Path(temp_dir) / "alignment.fasta"

        if self.alignment_tool == "mafft":
            cmd = ["mafft", "--auto", str(input_file)]
        elif self.alignment_tool == "muscle":
            cmd = ["muscle", "-in", str(input_file), "-out", str(output_file)]
        else:
            raise ValueError(f"Unsupported alignment tool: {self.alignment_tool}")

        try:
            if self.alignment_tool == "mafft":
                result = subprocess.run(cmd, capture_output=True, text=True, check=True)
                with open(output_file, "w") as f:
                    f.write(result.stdout)
            else:
                subprocess.run(cmd, check=True)

            return str(output_file)

        except subprocess.CalledProcessError as e:
            print(f"Alignment failed: {e}")
            return None
        except FileNotFoundError:
            print(
                f"Alignment tool '{self.alignment_tool}' not found. Please install it."
            )
            return None

    def get_reference_region(self, region: str) -> str:
        """Get reference sequence for specified HVS region."""
        if region not in self.hvs_regions:
            raise ValueError(f"Unknown region: {region}")

        start = self.hvs_regions[region]["start"] - 1  # Convert to 0-based
        end = self.hvs_regions[region]["end"]

        return self.reference_seq[start:end]

    def compare_sequences_direct(self, sample_seq: str, region: str) -> List[str]:
        """
        Direct comparison method (original).

        Args:
            sample_seq: Sample sequence
            region: HVS region

        Returns:
            List of mutations in HSD format
        """
        reference_region = self.get_reference_region(region)
        mutations = []
        start_pos = self.hvs_regions[region]["start"]

        # Simple position-by-position comparison
        min_len = min(len(sample_seq), len(reference_region))

        for i in range(min_len):
            if sample_seq[i] != reference_region[i]:
                pos = start_pos + i
                ref_base = reference_region[i]
                sample_base = sample_seq[i]
                mutations.append(f"{pos}{ref_base}>{sample_base}")

        # Handle length differences
        if len(sample_seq) > len(reference_region):
            # Insertion
            extra_bases = sample_seq[len(reference_region) :]
            pos = start_pos + len(reference_region)
            mutations.append(f"{pos}ins{extra_bases}")
        elif len(sample_seq) < len(reference_region):
            # Deletion
            deleted_bases = reference_region[len(sample_seq) :]
            pos = start_pos + len(sample_seq)
            mutations.append(f"{pos}del{deleted_bases}")

        return mutations

    def compare_sequences_aligned(self, sample_seq: str, region: str) -> List[str]:
        """
        Reference-aligned comparison method (enhanced).

        Args:
            sample_seq: Sample sequence
            region: HVS region

        Returns:
            List of mutations in HSD format
        """
        reference_region = self.get_reference_region(region)

        # Create temporary directory for alignment
        with tempfile.TemporaryDirectory() as temp_dir:
            # Align sample sequence to reference
            alignment_file = self.create_alignment(
                [reference_region, sample_seq], temp_dir
            )

            if not alignment_file:
                print("Warning: Alignment failed, falling back to direct comparison")
                return self.compare_sequences_direct(sample_seq, region)

            # Read alignment
            try:
                alignment = AlignIO.read(alignment_file, "fasta")
                ref_aligned = str(alignment[0].seq)
                sample_aligned = str(alignment[1].seq)

                return self.extract_mutations_from_alignment(
                    ref_aligned, sample_aligned, region
                )

            except Exception as e:
                print(f"Error reading alignment: {e}")
                return self.compare_sequences_direct(sample_seq, region)

    def extract_mutations_from_alignment(
        self, ref_aligned: str, sample_aligned: str, region: str
    ) -> List[str]:
        """Extract mutations from aligned sequences."""
        mutations = []
        start_pos = self.hvs_regions[region]["start"]
        ref_pos = start_pos - 1  # 0-based position in reference

        i = 0
        while i < len(ref_aligned) and i < len(sample_aligned):
            ref_base = ref_aligned[i]
            sample_base = sample_aligned[i]

            if ref_base == "-":
                # Insertion in sample
                insertion = ""
                j = i
                while j < len(sample_aligned) and ref_aligned[j] == "-":
                    insertion += sample_aligned[j]
                    j += 1
                if insertion:
                    mutations.append(f"{ref_pos + 1}ins{insertion}")
                i = j
                continue

            elif sample_base == "-":
                # Deletion in sample
                deletion = ""
                j = i
                while j < len(ref_aligned) and sample_aligned[j] == "-":
                    deletion += ref_aligned[j]
                    ref_pos += 1
                    j += 1
                if deletion:
                    mutations.append(f"{ref_pos}del{deletion}")
                i = j
                continue

            else:
                # Potential substitution
                ref_pos += 1
                if ref_base.upper() != sample_base.upper():
                    mutations.append(
                        f"{ref_pos}{ref_base.upper()}>{sample_base.upper()}"
                    )
                i += 1

        return mutations

    def convert_to_hsd(
        self, consensus_dir: str, output_file: str, method: str = "direct"
    ):
        """
        Convert consensus sequences to HSD format.

        Args:
            consensus_dir: Directory containing consensus FASTA files
            output_file: Output HSD file path
            method: Conversion method ("direct" or "aligned")
        """
        print(f"Converting consensus sequences to HSD format using {method} method...")

        # Load consensus sequences
        sequences = self.load_consensus_sequences(consensus_dir)

        if not sequences:
            print("No valid sequences found!")
            return

        # Process each sample
        with open(output_file, "w") as hsd_file:
            hsd_file.write("Sample\tMutations\n")

            for sample_id, regions in sequences.items():
                print(f"Processing sample: {sample_id}")
                all_mutations = []

                for region, sequence in regions.items():
                    print(f"  Processing {region}...")

                    if method == "aligned":
                        mutations = self.compare_sequences_aligned(sequence, region)
                    else:
                        mutations = self.compare_sequences_direct(sequence, region)

                    all_mutations.extend(mutations)
                    print(f"    Found {len(mutations)} mutations")

                # Write to HSD file
                mutations_str = ";".join(all_mutations) if all_mutations else "rCRS"
                hsd_file.write(f"{sample_id}\t{mutations_str}\n")
                print(f"  Total mutations for {sample_id}: {len(all_mutations)}")

        print(f"HSD conversion complete: {output_file}")
        print(f"Processed {len(sequences)} samples")


def main():
    """Main function for command-line usage."""
    if len(sys.argv) < 3:
        print(
            "Usage: python enhanced_hsd_converter.py <consensus_directory> <output_hsd_file> [--method=aligned]"
        )
        print("Methods:")
        print("  direct  - Direct position-by-position comparison (default)")
        print("  aligned - Reference-aligned comparison using external alignment tool")
        sys.exit(1)

    consensus_dir = sys.argv[1]
    output_file = sys.argv[2]

    # Parse method argument
    method = "direct"
    for arg in sys.argv[3:]:
        if arg.startswith("--method="):
            method = arg.split("=")[1]

    # Validate method
    if method not in ["direct", "aligned"]:
        print(f"Error: Unknown method '{method}'. Use 'direct' or 'aligned'.")
        sys.exit(1)

    # Create converter and run
    try:
        converter = EnhancedHSDConverter()
        converter.convert_to_hsd(consensus_dir, output_file, method)
    except Exception as e:
        print(f"Error: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
