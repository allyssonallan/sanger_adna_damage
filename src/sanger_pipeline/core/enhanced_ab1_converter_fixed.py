"""
Enhanced AB1 file converter with primer removal and trimming capabilities.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Tuple, Optional, Dict, Any
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import matplotlib.pyplot as plt

from ..utils.constants import DEFAULT_MIN_QUALITY, DEFAULT_MIN_SEQUENCE_LENGTH
from ..utils.helpers import validate_file_exists

logger = logging.getLogger(__name__)


class EnhancedAB1Converter:
    """
    Enhanced converter for AB1 files to FASTA format with primer removal, trimming, and quality filtering.
    """

    def __init__(
        self,
        min_quality: int = DEFAULT_MIN_QUALITY,
        min_sequence_length: int = DEFAULT_MIN_SEQUENCE_LENGTH,
        enable_primer_removal: bool = True,
        enable_quality_trimming: bool = True,
        quality_window_size: int = 10,
        quality_threshold_fraction: float = 0.75,
        adna_damage_mode: bool = True,
        adaptive_quality_threshold: bool = True,
        extremity_analysis_length: int = 30,
        custom_primers_forward: Optional[Dict[str, str]] = None,
        custom_primers_reverse: Optional[Dict[str, str]] = None,
    ):
        """
        Initialize enhanced AB1 converter.

        Args:
            min_quality: Minimum Phred quality score for base calling
            min_sequence_length: Minimum sequence length after filtering
            enable_primer_removal: Whether to remove primer sequences
            enable_quality_trimming: Whether to trim low-quality ends
            quality_window_size: Window size for quality trimming
            quality_threshold_fraction: Fraction of bases in window that must meet quality threshold
            adna_damage_mode: Enable ancient DNA damage-specific optimizations
            adaptive_quality_threshold: Dynamically adjust quality thresholds based on sequence characteristics
            extremity_analysis_length: Length of sequence extremities to analyze for N abundance
            custom_primers_forward: Custom forward primers dict {region: sequence}
            custom_primers_reverse: Custom reverse primers dict {region: sequence}
        """
        self.min_quality = min_quality
        self.min_sequence_length = min_sequence_length
        self.enable_primer_removal = enable_primer_removal
        self.enable_quality_trimming = enable_quality_trimming
        self.quality_window_size = quality_window_size
        self.quality_threshold_fraction = quality_threshold_fraction
        self.adna_damage_mode = adna_damage_mode
        self.adaptive_quality_threshold = adaptive_quality_threshold
        self.extremity_analysis_length = extremity_analysis_length

        # Ancient DNA damage-specific parameters
        if self.adna_damage_mode:
            self.primer_similarity_threshold = 0.35  # Lower for degraded primers
            self.iupac_tolerance = True  # Handle ambiguous bases
            self.damage_position_bias = True  # Consider 5' and 3' damage patterns
            self.conservative_trimming = True  # More aggressive quality trimming
        else:
            self.primer_similarity_threshold = 0.8  # Higher for modern DNA
            self.iupac_tolerance = False
            self.damage_position_bias = False
            self.conservative_trimming = False

        # Default primer sequences for HVS regions (F/R pairs)
        self.default_primers = {
            "HVS1": {
                "forward": "CACCATTAGCACCCAAAGCT",
                "reverse": "TGATTTCACGGAGGATGGTG",  # Forward direction (will be converted to RC for matching)
            },
            "HVS2": {
                "forward": "GGTCTATCACCCTATTAACCAC",
                "reverse": "CTGTTAAAAGTGCATACCGCCA",  # Forward direction (will be converted to RC for matching)
            },
            "HVS3": {
                "forward": "CCGCTTCTGGCCACAGCACT",
                "reverse": "GGTGATGTGAGCCCGTCTAAAC",  # Forward direction (will be converted to RC for matching)
            },
        }

        # Set up primer pairs (combine default with custom if provided)
        self.primers = self._setup_primer_pairs(
            custom_primers_forward, custom_primers_reverse
        )

    def _reverse_complement(self, sequence: str) -> str:
        """
        Get reverse complement of DNA sequence.

        Args:
            sequence: DNA sequence string

        Returns:
            Reverse complement sequence
        """
        complement_map = {
            "A": "T",
            "T": "A",
            "C": "G",
            "G": "C",
            "R": "Y",
            "Y": "R",
            "S": "S",
            "W": "W",
            "K": "M",
            "M": "K",
            "B": "V",
            "V": "B",
            "D": "H",
            "H": "D",
            "N": "N",
        }
        return "".join(
            complement_map.get(base.upper(), base) for base in sequence[::-1]
        )

    def _setup_primer_pairs(
        self,
        custom_forward: Optional[Dict[str, str]] = None,
        custom_reverse: Optional[Dict[str, str]] = None,
    ) -> Dict[str, Dict[str, str]]:
        """
        Set up primer pairs with proper forward/reverse complement handling.

        Args:
            custom_forward: Custom forward primers {region: sequence}
            custom_reverse: Custom reverse primers {region: sequence}

        Returns:
            Primer dictionary with forward and reverse_complement for each region
        """
        primers = {}

        # Start with default primers
        for region, primer_pair in self.default_primers.items():
            if "forward" not in primer_pair or "reverse" not in primer_pair:
                raise ValueError(
                    f"Missing forward or reverse primer for region {region}: {primer_pair}"
                )
            primers[region] = {
                "forward": primer_pair["forward"],
                "reverse_complement": self._reverse_complement(primer_pair["reverse"]),
                "reverse_original": primer_pair[
                    "reverse"
                ],  # Keep original for reference
            }

        # Add custom forward primers
        if custom_forward:
            for region, forward_seq in custom_forward.items():
                if region not in primers:
                    primers[region] = {}
                primers[region]["forward"] = forward_seq.upper()

        # Add custom reverse primers
        if custom_reverse:
            for region, reverse_seq in custom_reverse.items():
                if region not in primers:
                    primers[region] = {}
                primers[region]["reverse_complement"] = self._reverse_complement(
                    reverse_seq.upper()
                )
                primers[region]["reverse_original"] = reverse_seq.upper()

        return primers

    def detect_primer_orientation(self, sequence: str, region: str) -> Dict[str, Any]:
        """
        Detect primer orientation and quality in sequence.

        Args:
            sequence: DNA sequence to analyze
            region: HVS region to check

        Returns:
            Dictionary with primer detection results
        """
        if region not in self.primers:
            return {"detected": False, "reason": "unknown_region"}

        seq_upper = self._normalize_sequence(sequence)
        primers = self.primers[region]

        forward_primer = primers.get("forward", "")
        reverse_complement = primers.get("reverse_complement", "")

        # Check forward primer at beginning
        forward_score = self._find_best_primer_match(
            seq_upper, forward_primer, from_start=True
        )
        forward_pos = self._find_primer_position(
            seq_upper, forward_primer, from_start=True
        )

        # Check reverse complement at end
        reverse_score = self._find_best_primer_match(
            seq_upper, reverse_complement, from_start=False
        )
        reverse_pos = self._find_primer_position(
            seq_upper, reverse_complement, from_start=False
        )

        # Determine orientation
        if forward_score >= 0.4 and reverse_score >= 0.4:
            orientation = "F_to_R"  # Forward to Reverse (normal)
        elif forward_score >= 0.4:
            orientation = "F_only"  # Forward primer only
        elif reverse_score >= 0.4:
            orientation = "R_only"  # Reverse primer only
        else:
            orientation = "unknown"

        return {
            "detected": forward_score >= 0.4 or reverse_score >= 0.4,
            "orientation": orientation,
            "forward_score": forward_score,
            "reverse_score": reverse_score,
            "forward_position": forward_pos,
            "reverse_position": reverse_pos,
            "forward_primer": forward_primer,
            "reverse_complement": reverse_complement,
            "region": region,
        }

    def _normalize_sequence(self, sequence: str) -> str:
        """
        Normalize sequence by converting IUPAC codes and converting to uppercase.

        Args:
            sequence: Input DNA sequence

        Returns:
            Normalized sequence
        """
        seq_upper = sequence.upper()

        # Convert IUPAC ambiguous codes to N for better matching (optional)
        if self.iupac_tolerance:
            iupac_map = {
                "R": "N",
                "Y": "N",
                "S": "N",
                "W": "N",
                "K": "N",
                "M": "N",
                "B": "N",
                "D": "N",
                "H": "N",
                "V": "N",
            }
            for iupac, replacement in iupac_map.items():
                seq_upper = seq_upper.replace(iupac, replacement)

        return seq_upper

    def _find_primer_position(
        self, sequence: str, primer: str, from_start: bool = True
    ) -> int:
        """
        Find the position of best primer match.

        Args:
            sequence: Target sequence
            primer: Primer sequence to find
            from_start: If True, search from start; if False, search from end

        Returns:
            Position of best match (-1 if not found)
        """
        best_score = 0.0
        best_position = -1
        search_range = min(len(primer) + 10, len(sequence))

        if from_start:
            # Search forward primer from the beginning
            for i in range(search_range):
                window = sequence[i : i + len(primer)]
                if len(window) == len(primer):
                    score = self._calculate_match_score(window, primer)
                    if score > best_score and score >= 0.4:
                        best_score = score
                        best_position = i
        else:
            # Search reverse primer from the end
            for i in range(search_range):
                start_pos = len(sequence) - len(primer) - i
                if start_pos >= 0:
                    window = sequence[start_pos : start_pos + len(primer)]
                    if len(window) == len(primer):
                        score = self._calculate_match_score(window, primer)
                        if score > best_score and score >= 0.4:
                            best_score = score
                            best_position = start_pos

        return best_position

    def remove_primers_with_orientation(
        self, sequence: str, region: str
    ) -> Tuple[str, Dict]:
        """
        Remove primers with orientation detection and detailed reporting.

        Args:
            sequence: Input sequence
            region: HVS region

        Returns:
            Tuple of (cleaned_sequence, removal_info)
        """
        if not self.enable_primer_removal or region not in self.primers:
            return sequence, {
                "primers_removed": False,
                "reason": "disabled_or_unknown_region",
            }

        # Detect primer orientation
        primer_info = self.detect_primer_orientation(sequence, region)

        if not primer_info["detected"]:
            return sequence, {
                "primers_removed": False,
                "reason": "no_primers_detected",
                "primer_info": primer_info,
            }

        cleaned_seq = self._normalize_sequence(sequence)
        removal_info = {
            "primers_removed": False,
            "forward_removed": False,
            "reverse_removed": False,
            "original_length": len(sequence),
            "primer_info": primer_info,
        }

        # Remove forward primer if detected
        if primer_info["forward_score"] >= 0.4 and primer_info["forward_position"] >= 0:
            forward_end = primer_info["forward_position"] + len(
                primer_info["forward_primer"]
            )
            cleaned_seq = cleaned_seq[forward_end:]
            removal_info["forward_removed"] = True
            removal_info["forward_removed_length"] = forward_end

        # Remove reverse primer if detected
        if primer_info["reverse_score"] >= 0.4 and primer_info["reverse_position"] >= 0:
            # Adjust position if forward was removed
            reverse_pos = primer_info["reverse_position"]
            if removal_info["forward_removed"]:
                reverse_pos -= removal_info.get("forward_removed_length", 0)

            if reverse_pos >= 0 and reverse_pos < len(cleaned_seq):
                cleaned_seq = cleaned_seq[:reverse_pos]
                removal_info["reverse_removed"] = True

        removal_info["primers_removed"] = (
            removal_info["forward_removed"] or removal_info["reverse_removed"]
        )
        removal_info["final_length"] = len(cleaned_seq)
        removal_info["length_reduction"] = (
            removal_info["original_length"] - removal_info["final_length"]
        )

        return cleaned_seq, removal_info

    def analyze_extremity_n_abundance(
        self, sequence: str, stage: str = "unknown"
    ) -> Dict:
        """
        Analyze N abundance at sequence extremities (5' and 3' ends).

        Args:
            sequence: DNA sequence to analyze
            stage: Processing stage for tracking

        Returns:
            Dictionary with N abundance statistics
        """
        seq_len = len(sequence)
        analysis_len = min(self.extremity_analysis_length, seq_len // 2)

        if seq_len < analysis_len * 2:
            return {
                "stage": stage,
                "sequence_length": seq_len,
                "analysis_length": analysis_len,
                "five_prime_n_count": 0,
                "three_prime_n_count": 0,
                "five_prime_n_fraction": 0.0,
                "three_prime_n_fraction": 0.0,
                "total_n_count": sequence.count("N"),
                "total_n_fraction": (
                    sequence.count("N") / seq_len if seq_len > 0 else 0.0
                ),
                "damage_pattern": "insufficient_length",
            }

        # Analyze 5' end (first 30bp)
        five_prime_seq = sequence[:analysis_len].upper()
        five_prime_n_count = five_prime_seq.count("N")
        five_prime_n_fraction = five_prime_n_count / analysis_len

        # Analyze 3' end (last 30bp)
        three_prime_seq = sequence[-analysis_len:].upper()
        three_prime_n_count = three_prime_seq.count("N")
        three_prime_n_fraction = three_prime_n_count / analysis_len

        # Overall sequence statistics
        total_n_count = sequence.upper().count("N")
        total_n_fraction = total_n_count / seq_len

        # Determine damage pattern
        damage_pattern = self._classify_damage_pattern(
            five_prime_n_fraction, three_prime_n_fraction, total_n_fraction
        )

        return {
            "stage": stage,
            "sequence_length": seq_len,
            "analysis_length": analysis_len,
            "five_prime_n_count": five_prime_n_count,
            "three_prime_n_count": three_prime_n_count,
            "five_prime_n_fraction": five_prime_n_fraction,
            "three_prime_n_fraction": three_prime_n_fraction,
            "five_prime_sequence": five_prime_seq,
            "three_prime_sequence": three_prime_seq,
            "total_n_count": total_n_count,
            "total_n_fraction": total_n_fraction,
            "damage_pattern": damage_pattern,
        }

    def _classify_damage_pattern(
        self, five_prime_n: float, three_prime_n: float, total_n: float
    ) -> str:
        """
        Classify DNA damage pattern based on N distribution.

        Args:
            five_prime_n: Fraction of N's at 5' end
            three_prime_n: Fraction of N's at 3' end
            total_n: Overall fraction of N's

        Returns:
            Damage pattern classification
        """
        # Ancient DNA typically shows 5' > 3' damage (C→T at 5', G→A at 3')
        if five_prime_n > 0.4 and three_prime_n > 0.4:
            return "severe_degradation"
        elif five_prime_n > 0.2 or three_prime_n > 0.2:
            if five_prime_n > three_prime_n * 1.5:
                return "five_prime_bias"  # Typical aDNA pattern
            elif three_prime_n > five_prime_n * 1.5:
                return "three_prime_bias"
            else:
                return "moderate_degradation"
        elif total_n > 0.1:
            return "scattered_damage"
        else:
            return "minimal_damage"

    def adjust_parameters_for_damage(self, n_analysis: Dict) -> Dict:
        """
        Dynamically adjust processing parameters based on damage patterns.

        Args:
            n_analysis: N abundance analysis results

        Returns:
            Dictionary of adjusted parameters
        """
        if not self.adaptive_quality_threshold:
            return {"adjusted": False}

        damage_pattern = n_analysis.get("damage_pattern", "minimal_damage")

        # Base parameters
        adjusted_params = {
            "adjusted": True,
            "original_min_quality": self.min_quality,
            "original_window_size": self.quality_window_size,
            "original_threshold_fraction": self.quality_threshold_fraction,
            "damage_pattern": damage_pattern,
        }

        # Adjust based on damage severity
        if damage_pattern == "severe_degradation":
            # Very lenient for severely damaged DNA
            adjusted_params["min_quality"] = max(15, self.min_quality - 5)
            adjusted_params["quality_window_size"] = max(
                5, self.quality_window_size - 3
            )
            adjusted_params["quality_threshold_fraction"] = max(
                0.5, self.quality_threshold_fraction - 0.2
            )
            adjusted_params["primer_similarity_threshold"] = 0.25

        elif damage_pattern in [
            "five_prime_bias",
            "three_prime_bias",
            "moderate_degradation",
        ]:
            # Moderately lenient for damaged DNA
            adjusted_params["min_quality"] = max(17, self.min_quality - 3)
            adjusted_params["quality_window_size"] = max(
                7, self.quality_window_size - 2
            )
            adjusted_params["quality_threshold_fraction"] = max(
                0.6, self.quality_threshold_fraction - 0.1
            )
            adjusted_params["primer_similarity_threshold"] = 0.35

        elif damage_pattern == "scattered_damage":
            # Slightly more lenient
            adjusted_params["min_quality"] = max(18, self.min_quality - 2)
            adjusted_params["quality_window_size"] = self.quality_window_size
            adjusted_params["quality_threshold_fraction"] = max(
                0.65, self.quality_threshold_fraction - 0.05
            )
            adjusted_params["primer_similarity_threshold"] = 0.4

        else:  # minimal_damage
            # Use standard parameters
            adjusted_params["min_quality"] = self.min_quality
            adjusted_params["quality_window_size"] = self.quality_window_size
            adjusted_params["quality_threshold_fraction"] = (
                self.quality_threshold_fraction
            )
            adjusted_params["primer_similarity_threshold"] = (
                0.8 if not self.adna_damage_mode else 0.5
            )

        return adjusted_params

    def _create_damage_summary(self, processing_stats: Dict) -> Dict:
        """
        Create a comprehensive damage assessment summary.

        Args:
            processing_stats: Processing statistics dictionary

        Returns:
            Damage assessment summary
        """
        original = processing_stats.get("n_abundance_original", {})
        after_primers = processing_stats.get("n_abundance_after_primers", {})
        after_trimming = processing_stats.get("n_abundance_after_trimming", {})
        final = processing_stats.get("n_abundance_final", {})

        return {
            "damage_progression": {
                "original_pattern": original.get("damage_pattern", "unknown"),
                "after_primers_pattern": after_primers.get("damage_pattern", "unknown"),
                "after_trimming_pattern": after_trimming.get(
                    "damage_pattern", "unknown"
                ),
                "final_pattern": final.get("damage_pattern", "unknown"),
            },
            "n_fraction_changes": {
                "original_total": original.get("total_n_fraction", 0.0),
                "after_primers_total": after_primers.get("total_n_fraction", 0.0),
                "after_trimming_total": after_trimming.get("total_n_fraction", 0.0),
                "final_total": final.get("total_n_fraction", 0.0),
            },
            "extremity_analysis": {
                "original_5prime": original.get("five_prime_n_fraction", 0.0),
                "original_3prime": original.get("three_prime_n_fraction", 0.0),
                "final_5prime": final.get("five_prime_n_fraction", 0.0),
                "final_3prime": final.get("three_prime_n_fraction", 0.0),
                "five_prime_improvement": original.get("five_prime_n_fraction", 0.0)
                - final.get("five_prime_n_fraction", 0.0),
                "three_prime_improvement": original.get("three_prime_n_fraction", 0.0)
                - final.get("three_prime_n_fraction", 0.0),
            },
            "processing_effectiveness": {
                "primer_removal_effective": after_primers.get("total_n_fraction", 1.0)
                < original.get("total_n_fraction", 1.0),
                "trimming_effective": after_trimming.get("total_n_fraction", 1.0)
                < after_primers.get("total_n_fraction", 1.0),
                "overall_improvement": original.get("total_n_fraction", 1.0)
                - final.get("total_n_fraction", 1.0),
            },
        }

    def process_ab1_file_enhanced(
        self,
        ab1_file: Path,
        fasta_output: Path,
        processed_output: Path,
        plot_output: Path,
    ) -> Tuple[SeqRecord, Optional[SeqRecord], Dict]:
        """
        Complete enhanced processing of an AB1 file with damage analysis.

        Args:
            ab1_file: Path to input AB1 file
            fasta_output: Path to raw FASTA output
            processed_output: Path to processed FASTA output
            plot_output: Path to quality plot output

        Returns:
            Tuple of (raw_record, processed_record or None if too short, processing_stats)
        """
        logger.info(f"Enhanced processing AB1 file: {ab1_file}")

        # Convert to FASTA
        raw_record = self.convert_to_fasta(ab1_file, fasta_output)

        if "phred_quality" not in raw_record.letter_annotations:
            logger.warning(f"No quality scores found in {ab1_file}")
            return raw_record, None, {"error": "No quality scores"}

        sequence = str(raw_record.seq)
        qualities = raw_record.letter_annotations["phred_quality"]

        processing_stats = {
            "original_length": len(sequence),
            "hvs_region": None,
            "primers_removed": False,
            "quality_trimmed": (0, 0),
            "final_length": 0,
            "valid_bases": 0,
            "n_abundance_original": {},
            "n_abundance_after_primers": {},
            "n_abundance_after_trimming": {},
            "n_abundance_final": {},
            "dynamic_parameters": {},
            "damage_assessment": {},
        }

        # Step 0: Analyze original sequence N abundance
        n_analysis_original = self.analyze_extremity_n_abundance(sequence, "original")
        processing_stats["n_abundance_original"] = n_analysis_original

        # Adjust parameters based on initial damage assessment
        if self.adaptive_quality_threshold:
            adjusted_params = self.adjust_parameters_for_damage(n_analysis_original)
            processing_stats["dynamic_parameters"] = adjusted_params

            # Apply adjusted parameters for this sequence
            if adjusted_params.get("adjusted", False):
                current_min_quality = adjusted_params.get(
                    "min_quality", self.min_quality
                )
                current_window_size = adjusted_params.get(
                    "quality_window_size", self.quality_window_size
                )
                current_threshold_fraction = adjusted_params.get(
                    "quality_threshold_fraction", self.quality_threshold_fraction
                )
                logger.info(
                    f"Adjusted parameters for {n_analysis_original['damage_pattern']}: Q{current_min_quality}, window={current_window_size}, fraction={current_threshold_fraction:.2f}"
                )
            else:
                current_min_quality = self.min_quality
                current_window_size = self.quality_window_size
                current_threshold_fraction = self.quality_threshold_fraction
        else:
            current_min_quality = self.min_quality
            current_window_size = self.quality_window_size
            current_threshold_fraction = self.quality_threshold_fraction

        # Step 1: Detect HVS region
        hvs_region = self.detect_hvs_region(sequence)
        processing_stats["hvs_region"] = hvs_region
        logger.info(f"Detected HVS region: {hvs_region}")

        # Step 2: Remove primers with orientation detection
        if hvs_region and self.enable_primer_removal:
            primer_sequence, primer_removal_info = self.remove_primers_with_orientation(
                sequence, hvs_region
            )
            processing_stats["primers_removed"] = primer_removal_info.get(
                "primers_removed", False
            )
            processing_stats["primer_removal_info"] = primer_removal_info

            if primer_removal_info.get("primers_removed", False):
                sequence = primer_sequence  # Update sequence for further processing
                logger.info(
                    f"Primers removed: F={primer_removal_info.get('forward_removed', False)}, R={primer_removal_info.get('reverse_removed', False)}"
                )
                logger.info(
                    f"Length reduction: {primer_removal_info.get('length_reduction', 0)} bases"
                )
        else:
            processing_stats["primers_removed"] = False
            processing_stats["primer_removal_info"] = {
                "primers_removed": False,
                "reason": "disabled_or_no_region",
            }

        # Analyze N abundance after primer removal
        n_analysis_after_primers = self.analyze_extremity_n_abundance(
            sequence, "after_primers"
        )
        processing_stats["n_abundance_after_primers"] = n_analysis_after_primers

        # Step 3: Quality filtering (replace low-quality bases with N)
        filtered_seq = "".join(
            [
                base if qual >= current_min_quality else "N"
                for base, qual in zip(sequence, qualities)
            ]
        )

        # Step 4: Validate final sequence length
        valid_bases = sum(1 for base in filtered_seq if base in "ATCG")
        processing_stats["final_length"] = len(filtered_seq)
        processing_stats["valid_bases"] = valid_bases

        # Final N abundance analysis
        n_analysis_final = self.analyze_extremity_n_abundance(filtered_seq, "final")
        processing_stats["n_abundance_final"] = n_analysis_final

        # Create damage assessment summary
        processing_stats["damage_assessment"] = self._create_damage_summary(
            processing_stats
        )

        if valid_bases < self.min_sequence_length:
            logger.warning(
                f"Sequence {raw_record.id} too short after processing: {valid_bases} valid bases (minimum: {self.min_sequence_length})"
            )
            self.generate_quality_plot(raw_record, plot_output)
            return raw_record, None, processing_stats

        # Create processed record
        processed_record = SeqRecord(
            Seq(filtered_seq), id=raw_record.id, description=raw_record.description
        )

        # Write processed FASTA
        processed_output.parent.mkdir(parents=True, exist_ok=True)
        SeqIO.write(processed_record, processed_output, "fasta")
        logger.info(
            f"Wrote processed FASTA: {processed_output} ({valid_bases} valid bases)"
        )

        # Generate quality plot
        self.generate_quality_plot(raw_record, plot_output)

        logger.info(f"Enhanced processing completed: {ab1_file}")
        return raw_record, processed_record, processing_stats

    def convert_to_fasta(self, ab1_file: Path, output_file: Path) -> SeqRecord:
        """Convert AB1 file to FASTA format."""
        validate_file_exists(ab1_file, "AB1 file")

        try:
            record = SeqIO.read(ab1_file, "abi")
            logger.info(f"Read AB1 file: {ab1_file}")
        except Exception as e:
            raise ValueError(f"Failed to parse AB1 file {ab1_file}: {e}")

        record.id = output_file.name
        record.description = ""
        output_file.parent.mkdir(parents=True, exist_ok=True)
        SeqIO.write(record, output_file, "fasta")
        logger.info(f"Wrote FASTA file: {output_file}")
        return record

    def detect_hvs_region(self, sequence: str) -> Optional[str]:
        """Detect which HVS region based on primer presence."""
        seq_upper = sequence.upper()

        # Convert IUPAC ambiguous codes to N for better matching
        iupac_map = {
            "R": "N",
            "Y": "N",
            "S": "N",
            "W": "N",
            "K": "N",
            "M": "N",
            "B": "N",
            "D": "N",
            "H": "N",
            "V": "N",
        }
        for iupac, replacement in iupac_map.items():
            seq_upper = seq_upper.replace(iupac, replacement)

        best_region = None
        best_score = 0

        for hvs_region, primers in self.primers.items():
            forward_primer = primers["forward"]
            reverse_primer = primers[
                "reverse_complement"
            ]  # Use reverse complement for matching

            # Score based on best primer matches
            forward_score = self._find_best_primer_match(
                seq_upper, forward_primer, from_start=True
            )
            reverse_score = self._find_best_primer_match(
                seq_upper, reverse_primer, from_start=False
            )

            # Combined score (average of both primers)
            combined_score = (forward_score + reverse_score) / 2

            if (
                combined_score > best_score and combined_score > 0.35
            ):  # 35% minimum similarity
                best_score = combined_score
                best_region = hvs_region

        return best_region

    def _find_best_primer_match(
        self, sequence: str, primer: str, from_start: bool = True
    ) -> float:
        """Find best primer match allowing for significant mismatches."""
        best_score = 0.0
        search_range = min(len(primer) + 10, len(sequence))

        if from_start:
            for i in range(search_range):
                window = sequence[i : i + len(primer)]
                if len(window) == len(primer):
                    score = self._calculate_match_score(window, primer)
                    best_score = max(best_score, score)
        else:
            for i in range(search_range):
                start_pos = len(sequence) - len(primer) - i
                if start_pos >= 0:
                    window = sequence[start_pos : start_pos + len(primer)]
                    if len(window) == len(primer):
                        score = self._calculate_match_score(window, primer)
                        best_score = max(best_score, score)

        return best_score

    def _calculate_match_score(self, seq1: str, seq2: str) -> float:
        """Calculate match score between two sequences, ignoring N's."""
        if len(seq1) != len(seq2):
            return 0.0

        matches = 0
        valid_positions = 0

        for a, b in zip(seq1, seq2):
            if a != "N" and b != "N":  # Ignore positions with N
                valid_positions += 1
                if a == b:
                    matches += 1

        if valid_positions == 0:
            return 0.0

        return matches / valid_positions

    def generate_quality_plot(
        self,
        record: SeqRecord,
        output_file: Path,
        figure_size: Tuple[int, int] = (12, 4),
    ) -> None:
        """Generate quality score plot for the sequence."""
        if "phred_quality" not in record.letter_annotations:
            logger.warning("Record does not contain phred_quality annotations")
            return

        qualities = record.letter_annotations["phred_quality"]
        output_file.parent.mkdir(parents=True, exist_ok=True)

        plt.figure(figsize=figure_size)
        plt.plot(qualities, marker=".", linestyle="-", markersize=1)
        plt.axhline(
            self.min_quality,
            color="red",
            linestyle="--",
            label=f"Quality threshold ({self.min_quality})",
        )
        plt.title(f"Quality scores for {record.id}")
        plt.xlabel("Base position")
        plt.ylabel("Phred Quality")
        plt.legend()
        plt.grid(True, alpha=0.3)
        plt.tight_layout()
        plt.savefig(output_file, dpi=300, bbox_inches="tight")
        plt.close()
        logger.info(f"Generated quality plot: {output_file}")
