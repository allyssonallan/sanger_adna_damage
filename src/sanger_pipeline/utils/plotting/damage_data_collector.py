"""
Data collection utilities for damage visualization.

This module handles the collection and parsing of damage analysis data
from JSON files and directories.
"""

import json
import logging
from pathlib import Path
from typing import Dict, List, Any

logger = logging.getLogger(__name__)


class DamageDataCollector:
    """Collects damage analysis data from pipeline output."""

    def __init__(self, output_dir: Path):
        """
        Initialize damage data collector.

        Args:
            output_dir: Pipeline output directory
        """
        self.output_dir = Path(output_dir)
        self.damage_dir = self.output_dir / "damage_analysis"

    def collect_damage_data(self) -> List[Dict[str, Any]]:
        """
        Collect damage data from JSON files.

        Returns:
            List of damage data dictionaries
        """
        damage_data = []

        if not self.damage_dir.exists():
            logger.warning(f"Damage analysis directory not found: {self.damage_dir}")
            return damage_data

        json_files = list(self.damage_dir.glob("*.json"))

        for json_file in json_files:
            try:
                with open(json_file, "r") as f:
                    data = json.load(f)

                # Extract sample name from filename
                sample_name = json_file.stem.replace("_damage_results", "")

                # Parse damage data
                parsed_data = self._parse_damage_json(data, sample_name)
                if parsed_data:
                    damage_data.append(parsed_data)

            except Exception as e:
                logger.warning(f"Error reading damage file {json_file}: {e}")

        logger.info(f"Collected damage data for {len(damage_data)} samples")
        return damage_data

    def _parse_damage_json(self, data: Dict, sample_name: str) -> Dict[str, Any]:
        """
        Parse damage JSON data into standardized format.

        Args:
            data: Raw JSON data
            sample_name: Sample name

        Returns:
            Parsed damage data dictionary
        """
        try:
            damage_patterns = data.get("damage_patterns", {})
            sequence_quality = damage_patterns.get("sequence_quality", {})
            statistical_tests = damage_patterns.get("statistical_tests", {})

            return {
                "sample_name": sample_name,
                "damage_5_prime": damage_patterns.get("damage_5_prime", 0.0),
                "damage_3_prime": damage_patterns.get("damage_3_prime", 0.0),
                "overall_damage_rate": damage_patterns.get("overall_damage_rate", 0.0),
                "valid_percentage": sequence_quality.get("valid_percentage", 0.0),
                "n_percentage": sequence_quality.get("n_percentage", 0.0),
                "sequence_length": sequence_quality.get("sequence_length", 0),
                "p_value_5_prime": statistical_tests.get("p_value_5_prime", 1.0),
                "p_value_3_prime": statistical_tests.get("p_value_3_prime", 1.0),
                "confidence_5_prime": statistical_tests.get(
                    "confidence_5_prime", "low"
                ),
                "confidence_3_prime": statistical_tests.get(
                    "confidence_3_prime", "low"
                ),
                "position_data": damage_patterns.get("position_data", {}),
                "raw_data": data,  # Keep reference to full data
            }

        except Exception as e:
            logger.warning(f"Error parsing damage data for {sample_name}: {e}")
            return {}

    def get_dashboard_data(self, base_output_dir: str) -> Dict[str, Any]:
        """
        Generate dashboard data from damage analysis results.

        Args:
            base_output_dir: Base output directory path

        Returns:
            Dashboard data dictionary
        """
        base_path = Path(base_output_dir)
        damage_analysis_dir = base_path / "damage_analysis"

        if not damage_analysis_dir.exists():
            return {
                "has_data": False,
                "message": "No damage analysis directory found",
                "samples": [],
                "summary": {},
            }

        json_files = list(damage_analysis_dir.glob("*.json"))

        if not json_files:
            return {
                "has_data": False,
                "message": "No damage analysis JSON files found",
                "samples": [],
                "summary": {},
            }

        samples_data = []
        all_damage_5 = []
        all_damage_3 = []
        all_overall = []
        all_quality = []

        for json_file in json_files:
            try:
                with open(json_file, "r") as f:
                    data = json.load(f)

                sample_name = json_file.stem.replace("_damage_results", "")
                damage_patterns = data.get("damage_patterns", {})
                sequence_quality = damage_patterns.get("sequence_quality", {})
                statistical_tests = damage_patterns.get("statistical_tests", {})

                damage_5 = damage_patterns.get("damage_5_prime", 0.0)
                damage_3 = damage_patterns.get("damage_3_prime", 0.0)
                overall = damage_patterns.get("overall_damage_rate", 0.0)
                valid_pct = sequence_quality.get("valid_percentage", 0.0)

                sample_data = {
                    "name": sample_name,
                    "damage_5_prime": damage_5,
                    "damage_3_prime": damage_3,
                    "overall_damage": overall,
                    "quality_percentage": valid_pct,
                    "sequence_length": sequence_quality.get("sequence_length", 0),
                    "p_value_5": statistical_tests.get("p_value_5_prime", 1.0),
                    "p_value_3": statistical_tests.get("p_value_3_prime", 1.0),
                    "quality_tier": self._get_quality_tier(valid_pct),
                    "damage_tier": self._get_damage_tier(overall),
                    "confidence_level": self._get_confidence_level(
                        statistical_tests.get("p_value_5_prime", 1.0),
                        statistical_tests.get("p_value_3_prime", 1.0),
                    ),
                }

                samples_data.append(sample_data)
                all_damage_5.append(damage_5)
                all_damage_3.append(damage_3)
                all_overall.append(overall)
                all_quality.append(valid_pct)

            except Exception as e:
                logger.warning(f"Error processing {json_file}: {e}")

        # Calculate summary statistics
        summary = {}
        if samples_data:
            summary = {
                "total_samples": len(samples_data),
                "mean_damage_5": sum(all_damage_5) / len(all_damage_5),
                "mean_damage_3": sum(all_damage_3) / len(all_damage_3),
                "mean_overall_damage": sum(all_overall) / len(all_overall),
                "mean_quality": sum(all_quality) / len(all_quality),
                "high_damage_samples": len(
                    [s for s in samples_data if s["damage_tier"] in ["high", "extreme"]]
                ),
                "high_quality_samples": len(
                    [
                        s
                        for s in samples_data
                        if s["quality_tier"] in ["excellent", "good"]
                    ]
                ),
                "significant_samples": len(
                    [
                        s
                        for s in samples_data
                        if s["confidence_level"] in ["high", "very_high"]
                    ]
                ),
            }

        return {
            "has_data": True,
            "samples": samples_data,
            "summary": summary,
            "files_processed": len(json_files),
        }

    def _get_quality_tier(self, valid_percentage: float) -> str:
        """Classify sequence quality into tiers."""
        if valid_percentage >= 90:
            return "excellent"
        elif valid_percentage >= 75:
            return "good"
        elif valid_percentage >= 50:
            return "fair"
        else:
            return "poor"

    def _get_damage_tier(self, damage_rate: float) -> str:
        """Classify damage rate into tiers."""
        if damage_rate >= 30:
            return "extreme"
        elif damage_rate >= 20:
            return "high"
        elif damage_rate >= 10:
            return "moderate"
        elif damage_rate >= 5:
            return "low"
        else:
            return "minimal"

    def _get_confidence_level(self, p_value_5: float, p_value_3: float) -> str:
        """Determine statistical confidence level."""
        min_p = min(p_value_5, p_value_3)

        if min_p < 0.001:
            return "very_high"
        elif min_p < 0.01:
            return "high"
        elif min_p < 0.05:
            return "moderate"
        else:
            return "low"
