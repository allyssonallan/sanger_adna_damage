"""
Consensus sequence building from alignments.
"""

import logging
import subprocess
from pathlib import Path
from typing import List, Dict
from Bio import SeqIO, AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from ..utils.helpers import validate_file_exists


logger = logging.getLogger(__name__)


class ConsensusBuilder:
    """
    Build consensus sequences from aligned reads.
    """

    def __init__(self, alignment_tool: str = "mafft", alignment_params: str = "--auto"):
        """
        Initialize consensus builder.

        Args:
            alignment_tool: Alignment tool to use (default: mafft)
            alignment_params: Parameters for alignment tool
        """
        self.alignment_tool = alignment_tool
        self.alignment_params = alignment_params

    def reverse_complement_sequence(self, input_file: Path, output_file: Path) -> None:
        """
        Generate reverse complement of a sequence.

        Args:
            input_file: Input FASTA file
            output_file: Output FASTA file with reverse complement
        """
        validate_file_exists(input_file, "FASTA file")

        record = SeqIO.read(input_file, "fasta")
        record.seq = record.seq.reverse_complement()
        record.id += "_RC"
        record.description = "Reverse-complemented"

        # Ensure output directory exists
        output_file.parent.mkdir(parents=True, exist_ok=True)

        SeqIO.write(record, output_file, "fasta")
        logger.info(f"Generated reverse complement: {output_file}")

    def align_sequences(self, sequence_files: List[Path], output_file: Path) -> None:
        """
        Align multiple sequences using external alignment tool.

        Args:
            sequence_files: List of FASTA files to align
            output_file: Output alignment file
        """
        for seq_file in sequence_files:
            validate_file_exists(seq_file, "sequence file")

        # Create temporary combined file
        temp_file = output_file.parent / f"temp_{output_file.stem}.fasta"

        try:
            # Combine sequences
            with open(temp_file, "w") as temp_out:
                for seq_file in sequence_files:
                    with open(seq_file, "r") as seq_in:
                        temp_out.write(seq_in.read())

            # Ensure output directory exists
            output_file.parent.mkdir(parents=True, exist_ok=True)

            # Run alignment
            cmd = (
                [self.alignment_tool] + self.alignment_params.split() + [str(temp_file)]
            )

            with open(output_file, "w") as out_handle:
                subprocess.run(
                    cmd,
                    stdout=out_handle,
                    stderr=subprocess.PIPE,
                    text=True,
                    check=True,
                )

            logger.info(f"Aligned sequences: {output_file}")

        except subprocess.CalledProcessError as e:
            logger.error(f"Alignment failed: {e.stderr}")
            raise
        finally:
            # Clean up temporary file
            if temp_file.exists():
                temp_file.unlink()

    def build_consensus(
        self, alignment_file: Path, output_file: Path, consensus_id: str = "consensus"
    ) -> SeqRecord:
        """
        Build consensus sequence from alignment.

        Args:
            alignment_file: Input alignment file
            output_file: Output consensus FASTA file
            consensus_id: ID for consensus sequence

        Returns:
            Consensus SeqRecord
        """
        validate_file_exists(alignment_file, "alignment file")

        alignment = AlignIO.read(alignment_file, "fasta")

        consensus_seq = ""
        for i in range(alignment.get_alignment_length()):
            bases = [record.seq[i] for record in alignment]
            base_counts = self._count_bases(bases)
            consensus_base = self._determine_consensus_base(base_counts)
            consensus_seq += consensus_base

        # Create consensus record
        consensus_record = SeqRecord(
            Seq(consensus_seq), id=consensus_id, description="Consensus sequence"
        )

        # Ensure output directory exists
        output_file.parent.mkdir(parents=True, exist_ok=True)

        # Write consensus
        SeqIO.write(consensus_record, output_file, "fasta")
        logger.info(f"Built consensus: {output_file}")

        return consensus_record

    def _count_bases(self, bases: List[str]) -> Dict[str, int]:
        """
        Count occurrences of each base.

        Args:
            bases: List of bases at a position

        Returns:
            Dictionary with base counts
        """
        counts = {}
        for base in bases:
            if base not in ["N", "-"]:  # Ignore gaps and Ns
                counts[base] = counts.get(base, 0) + 1
        return counts

    def _determine_consensus_base(self, base_counts: Dict[str, int]) -> str:
        """
        Determine consensus base from counts.

        Args:
            base_counts: Dictionary with base counts

        Returns:
            Consensus base
        """
        if not base_counts:
            return "N"

        if len(base_counts) == 1:
            return list(base_counts.keys())[0]
        else:
            # If there's disagreement, return N
            return "N"

    def merge_regions(
        self, region_files: List[Path], output_file: Path, merged_id: str = "merged"
    ) -> SeqRecord:
        """
        Merge multiple sequence regions into a single sequence.

        Args:
            region_files: List of FASTA files containing regions to merge
            output_file: Output merged FASTA file
            merged_id: ID for merged sequence

        Returns:
            Merged SeqRecord
        """
        sequences = []

        for region_file in region_files:
            validate_file_exists(region_file, "region file")
            record = SeqIO.read(region_file, "fasta")
            sequences.append(str(record.seq))

        # Concatenate sequences
        merged_seq = "".join(sequences)

        # Create merged record
        merged_record = SeqRecord(
            Seq(merged_seq), id=merged_id, description="Merged sequence regions"
        )

        # Ensure output directory exists
        output_file.parent.mkdir(parents=True, exist_ok=True)

        # Write merged sequence
        SeqIO.write(merged_record, output_file, "fasta")
        logger.info(f"Merged regions: {output_file}")

        return merged_record

    def process_paired_reads(
        self,
        forward_file: Path,
        reverse_file: Path,
        alignment_output: Path,
        consensus_output: Path,
        sample_name: str,
    ) -> SeqRecord:
        """
        Process paired forward and reverse reads to generate consensus.

        Args:
            forward_file: Forward read FASTA file
            reverse_file: Reverse read FASTA file
            alignment_output: Output alignment file
            consensus_output: Output consensus file
            sample_name: Sample name for consensus

        Returns:
            Consensus SeqRecord
        """
        logger.info(f"Processing paired reads for sample: {sample_name}")

        # Create reverse complement of reverse read
        reverse_rc_file = reverse_file.parent / f"{reverse_file.stem}_rc.fasta"
        self.reverse_complement_sequence(reverse_file, reverse_rc_file)

        try:
            # Align sequences
            self.align_sequences([forward_file, reverse_rc_file], alignment_output)

            # Build consensus
            consensus_record = self.build_consensus(
                alignment_output,
                consensus_output,
                consensus_id=f"{sample_name}_consensus",
            )

            return consensus_record

        finally:
            # Clean up temporary reverse complement file
            if reverse_rc_file.exists():
                reverse_rc_file.unlink()
