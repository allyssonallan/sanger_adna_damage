#!/usr/bin/env python3
"""
HSD Diversity Analyzer

This tool analyzes HSD files to assess haplogroup diversity and detect
potential alignment artifacts or quality issues.

Author: Sanger aDNA Pipeline
"""

import logging
import argparse
from collections import Counter, defaultdict
import re

logger = logging.getLogger(__name__)


class HSDDiversityAnalyzer:
    """
    Analyzer for HSD files to assess diversity and quality.
    """

    def __init__(self):
        """Initialize the analyzer."""
        # Known haplogroup diagnostic positions (examples)
        self.diagnostic_positions = {
            "H": [16519],  # H haplogroup marker
            "U": [16189, 16270],  # U haplogroup markers
            "K": [16224, 16311],  # K haplogroup markers
            "T": [16126, 16294],  # T haplogroup markers
            "J": [16069, 16126],  # J haplogroup markers
            "V": [16298],  # V haplogroup marker
            "X": [16189, 16223],  # X haplogroup markers
        }

        # HVS regions
        self.hvs_regions = {
            "HVS1": (16024, 16365),
            "HVS2": (57, 372),
            "HVS3": (438, 574),
        }

    def parse_hsd_file(self, hsd_file: str) -> list:
        """
        Parse HSD file and extract sample data.

        Args:
            hsd_file: Path to HSD file

        Returns:
            List of sample dictionaries
        """
        samples = []

        with open(hsd_file, "r") as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue

                parts = line.split("\t")
                if len(parts) >= 4:
                    sample_id = parts[0]
                    regions = parts[1]
                    quality = parts[2]

                    # All parts from index 3 onwards are variants
                    variant_parts = parts[3:]

                    # Parse variants
                    variants = self._parse_variant_list(variant_parts)

                    sample_data = {
                        "id": sample_id,
                        "regions": regions,
                        "quality": quality,
                        "variants": variants,
                        "variant_count": len(variants),
                    }

                    samples.append(sample_data)

        return samples

    def _parse_variant_list(self, variant_parts: list) -> dict:
        """
        Parse list of variants into position-mutation dictionary.

        Args:
            variant_parts: List of variant strings (e.g., ["16026T", "16030A"])

        Returns:
            Dictionary of position -> mutation
        """
        variants = {}

        for variant in variant_parts:
            variant = variant.strip()
            if not variant or variant == "?":
                continue

            # Extract position and mutation
            match = re.match(r"(\d+)([ATGC])", variant)
            if match:
                position = int(match.group(1))
                mutation = match.group(2)
                variants[position] = mutation

        return variants

    def _parse_variants(self, variants_str: str) -> dict:
        """
        Parse variant string into position-mutation dictionary.

        Args:
            variants_str: String of variants (e.g., "16026T 16030A")

        Returns:
            Dictionary of position -> mutation
        """
        variants = {}

        if not variants_str or variants_str == "?":
            return variants

        # Split variants and parse each one
        variant_list = variants_str.split()

        for variant in variant_list:
            # Extract position and mutation
            match = re.match(r"(\d+)([ATGC])", variant)
            if match:
                position = int(match.group(1))
                mutation = match.group(2)
                variants[position] = mutation

        return variants

    def analyze_diversity(self, samples: list) -> dict:
        """
        Analyze diversity in the sample set.

        Args:
            samples: List of sample dictionaries

        Returns:
            Dictionary containing diversity analysis
        """
        analysis = {
            "total_samples": len(samples),
            "variant_statistics": {},
            "position_diversity": {},
            "sample_similarity": {},
            "quality_distribution": {},
            "potential_issues": [],
        }

        # Variant count statistics
        variant_counts = [s["variant_count"] for s in samples]
        analysis["variant_statistics"] = {
            "min_variants": min(variant_counts) if variant_counts else 0,
            "max_variants": max(variant_counts) if variant_counts else 0,
            "mean_variants": (
                sum(variant_counts) / len(variant_counts) if variant_counts else 0
            ),
            "total_unique_positions": len(self._get_all_positions(samples)),
        }

        # Position diversity analysis
        analysis["position_diversity"] = self._analyze_position_diversity(samples)

        # Sample similarity analysis
        analysis["sample_similarity"] = self._analyze_sample_similarity(samples)

        # Check for potential issues
        analysis["potential_issues"] = self._detect_potential_issues(samples)

        return analysis

    def _get_all_positions(self, samples: list) -> set:
        """Get all variant positions across all samples."""
        all_positions = set()
        for sample in samples:
            all_positions.update(sample["variants"].keys())
        return all_positions

    def _analyze_position_diversity(self, samples: list) -> dict:
        """
        Analyze diversity at each position.

        Args:
            samples: List of sample dictionaries

        Returns:
            Dictionary with position diversity analysis
        """
        position_mutations = defaultdict(list)

        # Collect all mutations at each position
        for sample in samples:
            for position, mutation in sample["variants"].items():
                position_mutations[position].append(mutation)

        diversity_stats = {}

        for position, mutations in position_mutations.items():
            mutation_counts = Counter(mutations)
            diversity_stats[position] = {
                "total_samples": len(mutations),
                "unique_mutations": len(mutation_counts),
                "mutations": dict(mutation_counts),
                "diversity_score": (
                    len(mutation_counts) / len(mutations) if mutations else 0
                ),
            }

        return diversity_stats

    def _analyze_sample_similarity(self, samples: list) -> dict:
        """
        Analyze similarity between samples.

        Args:
            samples: List of sample dictionaries

        Returns:
            Dictionary with similarity analysis
        """
        if len(samples) < 2:
            return {"message": "Need at least 2 samples for similarity analysis"}

        similarities = []

        for i in range(len(samples)):
            for j in range(i + 1, len(samples)):
                sample1 = samples[i]
                sample2 = samples[j]

                # Calculate Jaccard similarity
                variants1 = set(
                    f"{pos}{mut}" for pos, mut in sample1["variants"].items()
                )
                variants2 = set(
                    f"{pos}{mut}" for pos, mut in sample2["variants"].items()
                )

                intersection = len(variants1 & variants2)
                union = len(variants1 | variants2)

                jaccard_similarity = intersection / union if union > 0 else 0

                similarities.append(
                    {
                        "sample1": sample1["id"],
                        "sample2": sample2["id"],
                        "jaccard_similarity": jaccard_similarity,
                        "shared_variants": intersection,
                        "total_variants": union,
                    }
                )

        # Calculate overall statistics
        jaccard_scores = [s["jaccard_similarity"] for s in similarities]

        return {
            "pairwise_similarities": similarities,
            "mean_similarity": (
                sum(jaccard_scores) / len(jaccard_scores) if jaccard_scores else 0
            ),
            "min_similarity": min(jaccard_scores) if jaccard_scores else 0,
            "max_similarity": max(jaccard_scores) if jaccard_scores else 0,
        }

    def _detect_potential_issues(self, samples: list) -> list:
        """
        Detect potential quality issues or artifacts.

        Args:
            samples: List of sample dictionaries

        Returns:
            List of potential issues
        """
        issues = []

        # Check for very similar samples (potential contamination)
        similarities = self._analyze_sample_similarity(samples)
        if "pairwise_similarities" in similarities:
            high_similarity = [
                s
                for s in similarities["pairwise_similarities"]
                if s["jaccard_similarity"] > 0.8
            ]
            if high_similarity:
                issues.append(
                    {
                        "type": "High similarity",
                        "description": f"Found {len(high_similarity)} sample pairs with >80% similarity",
                        "details": high_similarity,
                    }
                )

        # Check for suspiciously low diversity
        if similarities.get("mean_similarity", 0) > 0.5:
            issues.append(
                {
                    "type": "Low diversity",
                    "description": f"Mean similarity is {similarities.get('mean_similarity', 0):.3f} (>0.5)",
                    "suggestion": "Check for contamination or alignment artifacts",
                }
            )

        # Check for samples with too many or too few variants
        variant_counts = [s["variant_count"] for s in samples]
        mean_variants = (
            sum(variant_counts) / len(variant_counts) if variant_counts else 0
        )

        for sample in samples:
            if sample["variant_count"] > mean_variants * 2:
                issues.append(
                    {
                        "type": "High variant count",
                        "sample": sample["id"],
                        "variant_count": sample["variant_count"],
                        "description": f"Sample has {sample['variant_count']} variants (mean: {mean_variants:.1f})",
                    }
                )
            elif sample["variant_count"] < mean_variants * 0.5:
                issues.append(
                    {
                        "type": "Low variant count",
                        "sample": sample["id"],
                        "variant_count": sample["variant_count"],
                        "description": f"Sample has {sample['variant_count']} variants (mean: {mean_variants:.1f})",
                    }
                )

        return issues

    def generate_report(self, hsd_file: str, output_file: str | None = None) -> str:
        """
        Generate comprehensive diversity report.

        Args:
            hsd_file: Path to HSD file
            output_file: Optional output file for report

        Returns:
            Report text
        """
        logger.info(f"Analyzing HSD file: {hsd_file}")

        # Parse HSD file
        samples = self.parse_hsd_file(hsd_file)

        # Analyze diversity
        analysis = self.analyze_diversity(samples)

        # Generate report
        report_lines = []
        report_lines.append("=" * 80)
        report_lines.append("HSD DIVERSITY ANALYSIS REPORT")
        report_lines.append("=" * 80)
        report_lines.append("")

        # Basic statistics
        report_lines.append("BASIC STATISTICS:")
        report_lines.append(f"  Total samples: {analysis['total_samples']}")
        report_lines.append("  Variants per sample:")
        report_lines.append(
            f"    Min: {analysis['variant_statistics']['min_variants']}"
        )
        report_lines.append(
            f"    Max: {analysis['variant_statistics']['max_variants']}"
        )
        report_lines.append(
            f"    Mean: {analysis['variant_statistics']['mean_variants']:.1f}"
        )
        report_lines.append(
            f"  Total unique positions: {analysis['variant_statistics']['total_unique_positions']}"
        )
        report_lines.append("")

        # Sample details
        report_lines.append("SAMPLE DETAILS:")
        for sample in samples:
            report_lines.append(f"  {sample['id']}: {sample['variant_count']} variants")
        report_lines.append("")

        # Diversity analysis
        report_lines.append("DIVERSITY ANALYSIS:")
        similarities = analysis["sample_similarity"]
        if "mean_similarity" in similarities:
            report_lines.append(
                f"  Mean sample similarity: {similarities['mean_similarity']:.3f}"
            )
            report_lines.append(
                f"  Min sample similarity: {similarities['min_similarity']:.3f}"
            )
            report_lines.append(
                f"  Max sample similarity: {similarities['max_similarity']:.3f}"
            )
        report_lines.append("")

        # Most variable positions
        position_div = analysis["position_diversity"]
        if position_div:
            sorted_positions = sorted(
                position_div.items(),
                key=lambda x: x[1]["diversity_score"],
                reverse=True,
            )
            report_lines.append("MOST VARIABLE POSITIONS (top 10):")
            for position, data in sorted_positions[:10]:
                mutations = ", ".join(
                    f"{mut}({count})" for mut, count in data["mutations"].items()
                )
                report_lines.append(
                    f"  Position {position}: {mutations} (diversity: {data['diversity_score']:.3f})"
                )
            report_lines.append("")

        # Potential issues
        issues = analysis["potential_issues"]
        if issues:
            report_lines.append("POTENTIAL ISSUES:")
            for issue in issues:
                report_lines.append(f"  - {issue['type']}: {issue['description']}")
                if "suggestion" in issue:
                    report_lines.append(f"    Suggestion: {issue['suggestion']}")
            report_lines.append("")
        else:
            report_lines.append("POTENTIAL ISSUES: None detected")
            report_lines.append("")

        # Recommendations
        report_lines.append("RECOMMENDATIONS:")
        if similarities.get("mean_similarity", 0) > 0.5:
            report_lines.append(
                "  - High similarity detected - check for contamination or alignment artifacts"
            )
            report_lines.append("  - Consider using stricter quality filters")
            report_lines.append("  - Verify reference sequence alignment")
        else:
            report_lines.append(
                "  - Good diversity detected - samples appear genetically distinct"
            )
            report_lines.append(
                "  - Consider uploading to HaploGrep for haplogroup classification"
            )

        report_lines.append("")
        report_lines.append("=" * 80)

        report_text = "\n".join(report_lines)

        # Write to file if specified
        if output_file:
            with open(output_file, "w") as f:
                f.write(report_text)
            logger.info(f"Report written to: {output_file}")

        return report_text


def main():
    """Main function for command line usage."""
    parser = argparse.ArgumentParser(description="HSD Diversity Analyzer")
    parser.add_argument("-i", "--input", required=True, help="Input HSD file")
    parser.add_argument("-o", "--output", help="Output report file")
    parser.add_argument("-v", "--verbose", action="store_true", help="Verbose logging")

    args = parser.parse_args()

    # Set up logging
    log_level = logging.DEBUG if args.verbose else logging.INFO
    logging.basicConfig(
        level=log_level, format="%(asctime)s - %(levelname)s - %(message)s"
    )

    # Create analyzer and generate report
    analyzer = HSDDiversityAnalyzer()
    report = analyzer.generate_report(args.input, args.output)

    # Print report to console
    print(report)


if __name__ == "__main__":
    main()
