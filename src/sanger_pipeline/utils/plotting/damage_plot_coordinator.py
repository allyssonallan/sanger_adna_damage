"""
Damage Plot Coordinator - Main class for orchestrating damage visualization.

This module coordinates between data collection and plot generation components
to create comprehensive damage analysis visualizations.
"""

import os
import argparse
from pathlib import Path
from typing import Dict, Any, Optional

from .damage_data_collector import DamageDataCollector
from .plot_generator import PlotGenerator


class DamagePlotCoordinator:
    """
    Main coordinator class for damage visualization.

    Orchestrates data collection and plot generation to create
    comprehensive damage analysis visualizations.
    """

    def __init__(self, output_dir: str):
        """Initialize coordinator with output directory."""
        self.output_dir = Path(output_dir)
        self.data_collector = DamageDataCollector(self.output_dir)
        self.plot_generator = PlotGenerator(self.output_dir / "plots")

    def generate_all_plots(self, output_dir: Optional[str] = None) -> Dict[str, str]:
        """
        Generate all damage analysis plots.

        Args:
            output_dir: Optional custom output directory

        Returns:
            Dictionary mapping plot types to their base64 encoded images
        """
        if output_dir is not None:
            self.output_dir = Path(output_dir)
            self.data_collector = DamageDataCollector(self.output_dir)
            self.plot_generator = PlotGenerator(self.output_dir / "plots")

        # Collect data
        damage_data = self.data_collector.collect_damage_data()

        # Generate plots
        plots = {}

        try:
            plots["damage_distribution"] = (
                self.plot_generator.create_damage_distribution_plot(damage_data)
            )
        except Exception as e:
            print(f"Error generating damage distribution plot: {e}")

        try:
            plots["correlation"] = self.plot_generator.create_damage_correlation_plot(
                damage_data
            )
        except Exception as e:
            print(f"Error generating correlation plot: {e}")

        try:
            plots["status_summary"] = self.plot_generator.create_status_summary_plot(
                damage_data
            )
        except Exception as e:
            print(f"Error generating status summary plot: {e}")

        try:
            plots["quality_damage"] = self.plot_generator.create_quality_damage_plot(
                damage_data
            )
        except Exception as e:
            print(f"Error generating quality damage plot: {e}")

        return plots

    def generate_dashboard_plots(self, base_output_dir: str) -> Dict[str, Any]:
        """
        Generate plots specifically for dashboard integration.

        Args:
            base_output_dir: Base output directory for analysis

        Returns:
            Dashboard-formatted plot data
        """
        dashboard_data = self.data_collector.get_dashboard_data(base_output_dir)

        return {
            "plots": self.generate_all_plots(),
            "data": dashboard_data,
            "summary": self._generate_summary_stats(dashboard_data),
        }

    def _generate_summary_stats(self, dashboard_data: Dict[str, Any]) -> Dict[str, Any]:
        """Generate summary statistics from dashboard data."""
        samples = dashboard_data.get("samples", [])

        if not samples:
            return {}

        # Calculate overall statistics
        total_samples = len(samples)

        # Quality statistics
        quality_counts = {}
        damage_counts = {}

        for sample in samples:
            quality_tier = sample.get("quality_tier", "unknown")
            damage_tier = sample.get("damage_tier", "unknown")

            quality_counts[quality_tier] = quality_counts.get(quality_tier, 0) + 1
            damage_counts[damage_tier] = damage_counts.get(damage_tier, 0) + 1

        # Average metrics
        valid_percentages = [s.get("valid_percentage", 0) for s in samples]
        damage_rates = [s.get("overall_damage", 0) for s in samples]

        avg_valid = (
            sum(valid_percentages) / len(valid_percentages) if valid_percentages else 0
        )
        avg_damage = sum(damage_rates) / len(damage_rates) if damage_rates else 0

        return {
            "total_samples": total_samples,
            "average_valid_percentage": round(avg_valid, 2),
            "average_damage_rate": round(avg_damage, 2),
            "quality_distribution": quality_counts,
            "damage_distribution": damage_counts,
        }


def main():
    """Command line interface for damage plot generation."""
    parser = argparse.ArgumentParser(description="Generate damage analysis plots")
    parser.add_argument("output_dir", help="Pipeline output directory")
    parser.add_argument(
        "--plots-dir", default="plots", help="Output directory for plots"
    )
    parser.add_argument(
        "--plot-type",
        choices=["all", "damage", "correlation", "status", "quality"],
        default="all",
        help="Type of plot to generate",
    )

    args = parser.parse_args()

    # Create plots directory if it doesn't exist
    plots_dir = Path(args.output_dir) / args.plots_dir
    plots_dir.mkdir(parents=True, exist_ok=True)

    # Initialize coordinator
    coordinator = DamagePlotCoordinator(args.output_dir)

    if args.plot_type == "all":
        plots = coordinator.generate_all_plots()
        print(f"Generated {len(plots)} plots in {plots_dir}")

        for plot_type, plot_data in plots.items():
            if plot_data:
                print(f"✓ {plot_type} plot generated successfully")
            else:
                print(f"✗ {plot_type} plot failed to generate")
    else:
        # Generate specific plot type
        data_collector = DamageDataCollector(Path(args.output_dir))
        plot_generator = PlotGenerator(plots_dir)

        damage_data = data_collector.collect_damage_data()

        plot_methods = {
            "damage": plot_generator.create_damage_distribution_plot,
            "correlation": plot_generator.create_damage_correlation_plot,
            "status": plot_generator.create_status_summary_plot,
            "quality": plot_generator.create_quality_damage_plot,
        }

        if args.plot_type in plot_methods:
            try:
                plot_data = plot_methods[args.plot_type](damage_data)
                if plot_data:
                    print(f"✓ {args.plot_type} plot generated successfully")
                else:
                    print(f"✗ {args.plot_type} plot failed to generate")
            except Exception as e:
                print(f"Error generating {args.plot_type} plot: {e}")


if __name__ == "__main__":
    main()
