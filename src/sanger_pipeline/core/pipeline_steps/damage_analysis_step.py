"""
Damage analysis step handler for the Sanger pipeline.

This module handles ancient DNA damage pattern analysis.
"""

import json
import logging
from pathlib import Path
from typing import Dict

from ..adna_damage_analyzer import ADNADamageAnalyzer

logger = logging.getLogger(__name__)


class DamageAnalysisStep:
    """Handles ancient DNA damage pattern analysis."""
    
    def __init__(self, damage_analyzer: ADNADamageAnalyzer, config: Dict):
        """
        Initialize the damage analysis step.
        
        Args:
            damage_analyzer: Configured ADNADamageAnalyzer instance
            config: Pipeline configuration dictionary
        """
        self.damage_analyzer = damage_analyzer
        self.config = config
    
    def execute(self, directories: Dict[str, Path]) -> Dict[str, int]:
        """
        Execute damage analysis step.
        
        Args:
            directories: Pipeline output directories
            
        Returns:
            Dictionary with analysis statistics
        """
        logger.info("Step 4: Analyzing aDNA damage patterns")

        damage_dir = directories["output"] / "damage_analysis"
        damage_dir.mkdir(exist_ok=True)

        # Get reference sequence
        ref_file = directories.get("ref", directories["output"].parent / "ref") / "rCRS.fasta"
        if not ref_file.exists():
            logger.warning(f"Reference file not found: {ref_file}. Skipping damage analysis.")
            return {"processed_samples": 0}

        processed_samples = 0
        consensus_files = list(directories["consensus"].glob("*_consensus.fasta"))

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
                self._save_damage_results(
                    damage_dir, sample_name, results, bootstrap_results, damage_assessment
                )

                logger.info(f"Damage analysis completed for {sample_name}")
                processed_samples += 1

            except Exception as e:
                logger.error(f"Failed to analyze damage patterns for {sample_name}: {e}")
                continue

        logger.info(f"Completed damage analysis for {processed_samples} samples")
        
        return {"processed_samples": processed_samples}
    
    def _save_damage_results(self, damage_dir: Path, sample_name: str, 
                           results, bootstrap_results, 
                           damage_assessment) -> None:
        """
        Save damage analysis results to JSON file.
        
        Args:
            damage_dir: Directory for damage analysis output
            sample_name: Name of the sample
            results: Damage pattern analysis results
            bootstrap_results: Bootstrap analysis results
            damage_assessment: Damage assessment results
        """
        results_file = damage_dir / f"{sample_name}_damage_results.json"
        
        with open(results_file, 'w') as f:
            json.dump({
                'damage_patterns': results,
                'bootstrap_analysis': bootstrap_results,
                'damage_assessment': damage_assessment
            }, f, indent=2)
