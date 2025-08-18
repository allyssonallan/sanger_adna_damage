"""
Refactored pipeline orchestrator using modular step handlers.

This module provides a clean, modular pipeline structure using separate step handlers.
"""

import logging
from pathlib import Path
from typing import Dict, Optional

from .ab1_converter import AB1Converter
from .consensus_builder import ConsensusBuilder
from .adna_damage_analyzer import ADNADamageAnalyzer
from ..utils.helpers import create_directories, load_config

from .pipeline_steps import (
    AB1ConversionStep, ConsensusStep, RegionMergingStep,
    DamageAnalysisStep, ReportGenerationStep, HSDConversionStep
)

logger = logging.getLogger(__name__)


class SangerPipelineRefactored:
    """
    Refactored main pipeline for processing Sanger sequencing data with aDNA damage analysis.
    
    Uses modular step handlers for improved maintainability and testability.
    """

    def __init__(
        self,
        input_dir: Path,
        output_dir: Path,
        config_file: Optional[Path] = None,
        min_quality: int = 20,
        min_sequence_length: int = 30,
        alignment: Optional[Dict] = None,
        alignment_tool: str = "mafft",
        alignment_params: str = "--auto",
    ):
        """
        Initialize the refactored Sanger pipeline.

        Args:
            input_dir: Directory containing AB1 files
            output_dir: Directory for output files
            config_file: Optional configuration file
            min_quality: Minimum Phred quality score
            min_sequence_length: Minimum sequence length after filtering
            alignment: Alignment configuration dict
            alignment_tool: Alignment tool name
            alignment_params: Alignment parameters
        """
        self.input_dir = Path(input_dir)
        self.output_dir = Path(output_dir)
        self.min_quality = min_quality
        self.min_sequence_length = min_sequence_length
        
        # Handle alignment configuration - prefer explicit parameters over dict
        if alignment:
            self.alignment = alignment
        else:
            self.alignment = {"tool": alignment_tool, "parameters": alignment_params}

        # Load configuration
        self.config = load_config(config_file) if config_file else {}

        # Override with config values if available
        quality_config = self.config.get('quality', {})
        self.min_quality = quality_config.get('min_phred_score', min_quality)
        self.min_sequence_length = quality_config.get('min_sequence_length', min_sequence_length)

        # Create output directories
        self.directories = create_directories(self.output_dir)

        # Initialize components
        self._initialize_components()
        
        # Initialize step handlers
        self._initialize_step_handlers()

        logger.info(f"Initialized refactored pipeline: {input_dir} -> {output_dir}")
        logger.info(f"Quality settings: min_quality={self.min_quality}, min_sequence_length={self.min_sequence_length}")

    def _initialize_components(self) -> None:
        """Initialize core pipeline components."""
        # Initialize components with sequence length filtering
        self.ab1_converter = AB1Converter(
            min_quality=self.min_quality, 
            min_sequence_length=self.min_sequence_length
        )
        self.consensus_builder = ConsensusBuilder()
        
        # Get damage threshold from config or use default
        damage_threshold = self.config.get('damage', {}).get('damage_threshold', 0.02)
        self.damage_analyzer = ADNADamageAnalyzer(min_damage_threshold=damage_threshold)

    def _initialize_step_handlers(self) -> None:
        """Initialize modular step handlers."""
        self.ab1_conversion_step = AB1ConversionStep(self.ab1_converter)
        self.consensus_step = ConsensusStep(self.consensus_builder)
        self.region_merging_step = RegionMergingStep()
        self.damage_analysis_step = DamageAnalysisStep(self.damage_analyzer, self.config)
        self.report_generation_step = ReportGenerationStep()
        self.hsd_conversion_step = HSDConversionStep(self.config)

    def run(self) -> Dict[str, Dict]:
        """
        Run the complete pipeline using modular step handlers.
        
        Returns:
            Dictionary containing results from each step
        """
        logger.info("Starting refactored Sanger sequencing pipeline")

        results = {}

        try:
            # Step 1: Convert AB1 files
            results['ab1_conversion'] = self.ab1_conversion_step.execute(
                self.input_dir, self.directories
            )

            # Step 2: Generate consensus sequences
            results['consensus_generation'] = self.consensus_step.execute(
                self.directories
            )

            # Step 3: Merge HVS regions
            results['region_merging'] = self.region_merging_step.execute(
                self.directories
            )

            # Step 4: Analyze aDNA damage patterns
            results['damage_analysis'] = self.damage_analysis_step.execute(
                self.directories
            )

            # Step 5: Generate comprehensive report
            results['report_generation'] = self.report_generation_step.execute(
                self.directories
            )

            # Step 6: Convert to HSD format
            results['hsd_conversion'] = self.hsd_conversion_step.execute(
                self.directories
            )

            logger.info("Refactored pipeline completed successfully")
            return results

        except Exception as e:
            logger.error(f"Pipeline failed: {e}")
            raise

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


    def _step_4_adna_damage_analysis(self, *args, **kwargs):
        """Backward compatibility method for old test."""
        return self.damage_analysis_step.execute(self.directories)


# Maintain backward compatibility
SangerPipeline = SangerPipelineRefactored