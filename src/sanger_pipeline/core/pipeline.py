"""
Main pipeline orchestrator for Sanger DNA damage analysis.
"""

import json
import logging
from pathlib import Path
from typing import Dict, Optional

from .ab1_converter import AB1Converter
from .consensus_builder import ConsensusBuilder
from .adna_damage_analyzer import ADNADamageAnalyzer
from ..utils.helpers import create_directories, load_config

logger = logging.getLogger(__name__)


class SangerPipeline:
    """
    Main pipeline for processing Sanger sequencing data with aDNA damage analysis.
    """

    def __init__(
        self,
        input_dir: Path,
        output_dir: Path,
        config_file: Optional[Path] = None,
        min_quality: int = 20,
        alignment: Optional[Dict] = None,
    ):
        """
        Initialize the Sanger pipeline.

        Args:
            input_dir: Directory containing AB1 files
            output_dir: Directory for output files
            config_file: Optional configuration file
            min_quality: Minimum Phred quality score
            alignment: Alignment configuration
        """
        self.input_dir = Path(input_dir)
        self.output_dir = Path(output_dir)
        self.min_quality = min_quality
        self.alignment = alignment or {"tool": "mafft", "parameters": "--auto"}

        # Load configuration
        self.config = load_config(config_file) if config_file else {}

        # Create output directories
        self.directories = create_directories(self.output_dir)

        # Initialize components
        self.ab1_converter = AB1Converter(min_quality=min_quality)
        self.consensus_builder = ConsensusBuilder()
        
        # Get damage threshold from config or use default
        damage_threshold = self.config.get('damage', {}).get('damage_threshold', 0.02)
        self.damage_analyzer = ADNADamageAnalyzer(min_damage_threshold=damage_threshold)

        logger.info(f"Initialized pipeline: {input_dir} -> {output_dir}")

    def run(self) -> None:
        """
        Run the complete pipeline.
        """
        logger.info("Starting Sanger sequencing pipeline")

        try:
            self._step_1_convert_ab1_files()
            self._step_2_align_and_consensus()
            self._step_3_merge_regions()
            self._step_4_adna_damage_analysis()
            self._step_5_generate_report()

            logger.info("Pipeline completed successfully")

        except Exception as e:
            logger.error(f"Pipeline failed: {e}")
            raise

    def _step_1_convert_ab1_files(self) -> None:
        """
        Step 1: Convert AB1 files to FASTA format with quality filtering.
        """
        logger.info("Step 1: Converting AB1 files to FASTA with quality filtering")

        ab1_files = list(self.input_dir.glob("*.ab1"))
        if not ab1_files:
            logger.warning(f"No AB1 files found in {self.input_dir}")
            return

        processed_files = 0
        for ab1_file in ab1_files:
            fasta_output = self.directories["fasta"] / f"{ab1_file.stem}.fasta"
            filtered_output = self.directories["filtered"] / f"{ab1_file.stem}_filtered.fasta"
            plot_output = self.directories["plots"] / f"{ab1_file.stem}_quality.png"

            try:
                # Convert AB1 to FASTA, filter by quality, and generate quality plot
                seq_record = self.ab1_converter.convert_to_fasta(ab1_file, fasta_output)
                filtered_record = self.ab1_converter.filter_by_quality(seq_record, filtered_output)
                self.ab1_converter.generate_quality_plot(seq_record, plot_output)
                processed_files += 1
            except Exception as e:
                logger.error(f"Failed to process {ab1_file}: {e}")
                continue

        logger.info(f"Converted {processed_files} AB1 files to FASTA and generated filtered versions")

    def _step_2_align_and_consensus(self) -> None:
        """
        Step 2: Align paired reads and generate consensus sequences for each HVS region.
        """
        logger.info("Step 2: Generating consensus sequences from paired reads by HVS region")

        fasta_files = list(self.directories["fasta"].glob("*.fasta"))
        if not fasta_files:
            logger.warning("No FASTA files found for consensus generation")
            return

        # Group files by sample and HVS region
        sample_hvs_groups = {}
        for fasta_file in fasta_files:
            name = fasta_file.stem
            
            # Parse filename pattern: SAMPLE_NAME_HVS#-Direction
            if "_HVS" in name and name.endswith(("-F", "-R")):
                # Split at HVS region
                hvs_start = name.rfind("_HVS")
                sample_name = name[:hvs_start]
                hvs_part = name[hvs_start+1:]  # Remove the '_'
                
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

        # Process each sample's HVS regions
        processed_samples = 0
        total_consensus = 0
        
        for sample_name, hvs_regions in sample_hvs_groups.items():
            logger.info(f"Processing sample {sample_name} with regions: {list(hvs_regions.keys())}")
            
            for hvs_region, files in hvs_regions.items():
                if "F" not in files or "R" not in files:
                    logger.warning(f"Missing paired reads for {sample_name} {hvs_region} (F:{files.get('F')} R:{files.get('R')})")
                    continue

                forward_file = files["F"]
                reverse_file = files["R"]

                # Output files named with sample and HVS region
                alignment_output = self.directories["aligned"] / f"{sample_name}_{hvs_region}_aligned.fasta"
                consensus_output = self.directories["consensus"] / f"{sample_name}_{hvs_region}_consensus.fasta"

                try:
                    self.consensus_builder.process_paired_reads(
                        forward_file, reverse_file, alignment_output, consensus_output, f"{sample_name}_{hvs_region}"
                    )
                    total_consensus += 1
                    logger.debug(f"Generated consensus for {sample_name} {hvs_region}")
                except Exception as e:
                    logger.error(f"Failed to process paired reads for {sample_name} {hvs_region}: {e}")
                    continue
            
            if hvs_regions:  # If this sample had any HVS regions processed
                processed_samples += 1

        logger.info(f"Generated consensus for {total_consensus} HVS regions across {processed_samples} samples")

    def _step_3_merge_regions(self) -> None:
        """
        Step 3: Merge available HVS regions into final consensus sequences.
        Creates combinations like HVS1_HVS2, HVS2_HVS3, HVS1_HVS2_HVS3 depending on availability.
        """
        logger.info("Step 3: Merging available HVS regions into final consensus sequences")

        consensus_files = list(self.directories["consensus"].glob("*_consensus.fasta"))
        if not consensus_files:
            logger.warning("No consensus files found for merging")
            return

        # Group consensus files by sample
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

        # Process each sample
        processed_samples = 0
        for sample_name, hvs_files in sample_consensus.items():
            available_regions = sorted(hvs_files.keys())
            logger.info(f"Sample {sample_name} has regions: {available_regions}")
            
            try:
                # Read all available HVS consensus sequences
                sequences = {}
                for hvs_region, consensus_file in hvs_files.items():
                    with open(consensus_file, 'r') as f:
                        # Read the sequence content (skip FASTA header)
                        lines = f.readlines()
                        sequence_lines = [line.strip() for line in lines if not line.startswith('>')]
                        sequences[hvs_region] = ''.join(sequence_lines)
                
                # Create merged sequence name based on available regions
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
                
                # Write merged consensus
                merged_output = self.directories["final"] / f"{merged_name}_merged.fasta"
                with open(merged_output, 'w') as f:
                    f.write(f">{merged_name}_consensus\n")
                    # Write sequence with line breaks every 80 characters
                    for i in range(0, len(merged_sequence), 80):
                        f.write(merged_sequence[i:i+80] + '\n')
                
                processed_samples += 1
                logger.debug(f"Created merged sequence for {sample_name}: {region_name}")
                
            except Exception as e:
                logger.error(f"Failed to merge regions for {sample_name}: {e}")
                continue

        logger.info(f"Created merged sequences for {processed_samples} samples")

    def _step_4_adna_damage_analysis(self) -> None:
        """
        Step 4: Ancient DNA damage pattern analysis.
        """
        logger.info("Step 4: Analyzing aDNA damage patterns")

        damage_dir = self.directories["output"] / "damage_analysis"
        damage_dir.mkdir(exist_ok=True)

        # Get reference sequence
        ref_file = self.directories.get("ref", self.output_dir.parent / "ref") / "rCRS.fasta"
        if not ref_file.exists():
            logger.warning(f"Reference file not found: {ref_file}. Skipping damage analysis.")
            return

        processed_samples = 0
        consensus_files = list(self.directories["consensus"].glob("*_consensus.fasta"))

        for consensus_file in consensus_files:
            sample_name = consensus_file.stem.replace("_consensus", "")
            
            try:
                # Analyze damage patterns
                results = self.damage_analyzer.analyze_sequence_damage(
                    consensus_file,
                    ref_file
                )

                # Generate damage plots using single file
                self.damage_analyzer.generate_damage_plots(
                    [consensus_file],
                    ref_file,
                    damage_dir
                )

                # Perform bootstrap analysis using single file
                bootstrap_results = self.damage_analyzer.bootstrap_damage_analysis(
                    [consensus_file],
                    ref_file,
                    iterations=self.config.get('bootstrap_iterations', 10000)
                )

                # Assess damage indicators
                damage_assessment = self.damage_analyzer.assess_damage_indicators(bootstrap_results)

                # Save results
                results_file = damage_dir / f"{sample_name}_damage_results.json"
                with open(results_file, 'w') as f:
                    json.dump({
                        'damage_patterns': results,
                        'bootstrap_analysis': bootstrap_results,
                        'damage_assessment': damage_assessment
                    }, f, indent=2)

                logger.info(f"Damage analysis completed for {sample_name}")
                processed_samples += 1

            except Exception as e:
                logger.error(f"Failed to analyze damage patterns for {sample_name}: {e}")
                continue

        logger.info(f"Completed damage analysis for {processed_samples} samples")

    def _step_5_generate_report(self) -> None:
        """
        Step 5: Generate comprehensive QC report.
        """
        logger.info("Step 5: Generating comprehensive QC report")

        # Include damage analysis results in the report
        damage_dir = self.directories["output"] / "damage_analysis"
        if damage_dir.exists():
            logger.info("Including aDNA damage analysis results in report")

        # Generate comprehensive QC report
        try:
            from ..utils.report_generator import QCReportGenerator
            
            report_generator = QCReportGenerator(self.directories["output"])
            report_file = report_generator.generate_report()
            
            logger.info(f"QC report generated successfully: {report_file}")
            
        except Exception as e:
            logger.error(f"Failed to generate QC report: {e}")
            logger.info("Report generation failed, but pipeline continues")

    def get_summary(self) -> Dict:
        """
        Get pipeline summary statistics.

        Returns:
            Dictionary with summary statistics
        """
        summary = {
            "input_directory": str(self.input_dir),
            "output_directory": str(self.output_dir),
            "ab1_files": len(list(self.input_dir.glob("*.ab1"))),
            "fasta_files": len(list(self.directories["fasta"].glob("*.fasta"))),
            "filtered_files": len(list(self.directories["filtered"].glob("*_filtered.fasta"))),
            "consensus_files": len(list(self.directories["consensus"].glob("*_consensus.fasta"))),
            "merged_files": len(list(self.directories["final"].glob("*_merged.fasta"))),
            "quality_plots": len(list(self.directories["plots"].glob("*_quality.png"))),
        }

        # Add damage analysis stats if available
        damage_dir = self.directories["output"] / "damage_analysis"
        if damage_dir.exists():
            summary["damage_results"] = len(list(damage_dir.glob("*_damage_results.json")))
            summary["damage_plots"] = len(list(damage_dir.glob("*_damage_plots.png")))

        return summary
