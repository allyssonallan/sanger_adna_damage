"""
Refactored damage plotting module using modular components.

This module provides a simplified interface to the modular damage plotting system,
maintaining backward compatibility while using the new component architecture.
"""

from pathlib import Path
from typing import Dict, Any, Optional

from .plotting.damage_plot_coordinator import DamagePlotCoordinator
from .plotting.damage_data_collector import DamageDataCollector
from .plotting.plot_generator import PlotGenerator


class DamagePlotGenerator:
    """Backward compatibility class for the old DamagePlotGenerator interface."""

    def __init__(self, output_dir: Optional[str] = None):
        """Initialize with output directory."""
        self.output_dir = output_dir
        self.coordinator = DamagePlotCoordinator(output_dir) if output_dir else None

    def generate_all_plots(self, output_dir: Optional[str] = None) -> Dict[str, str]:
        """Generate all plots with backward compatibility."""
        dir_to_use = output_dir or self.output_dir
        if not dir_to_use:
            raise ValueError("Output directory must be provided")

        coordinator = DamagePlotCoordinator(dir_to_use)
        return coordinator.generate_all_plots()

    def create_damage_plots(self, output_dir: str) -> Dict[str, str]:
        """Create damage plots (backward compatibility method)."""
        return self.generate_all_plots(output_dir)

    def get_dashboard_data(self, base_output_dir: str) -> Dict[str, Any]:
        """Get dashboard data for damage visualization."""
        coordinator = DamagePlotCoordinator(base_output_dir)
        return coordinator.generate_dashboard_plots(base_output_dir)

    def generate_comprehensive_damage_plots(self) -> Dict[str, str]:
        """Generate comprehensive damage plots for report embedding."""
        if not self.output_dir:
            return {}
        
        coordinator = DamagePlotCoordinator(self.output_dir)
        return coordinator.generate_all_plots()

    def generate_individual_sample_plots(self, output_dir: str) -> Dict[str, str]:
        """Generate individual sample plots."""
        coordinator = DamagePlotCoordinator(output_dir)
        return coordinator.generate_all_plots()


def create_damage_plots(output_dir: str) -> Dict[str, str]:
    """
    Create all damage analysis plots for the given output directory.

    Args:
        output_dir: Pipeline output directory containing damage analysis data

    Returns:
        Dictionary mapping plot types to their base64 encoded images
    """
    coordinator = DamagePlotCoordinator(output_dir)
    return coordinator.generate_all_plots()


def create_dashboard_data(base_output_dir: str) -> Dict[str, Any]:
    """
    Generate dashboard data for damage visualization.

    Args:
        base_output_dir: Base output directory for analysis

    Returns:
        Dashboard-formatted data including plots and statistics
    """
    coordinator = DamagePlotCoordinator(base_output_dir)
    return coordinator.generate_dashboard_plots(base_output_dir)


def create_damage_distribution_plot(output_dir: str) -> Optional[str]:
    """
    Create a damage distribution plot.

    Args:
        output_dir: Pipeline output directory

    Returns:
        Base64 encoded plot image
    """
    data_collector = DamageDataCollector(Path(output_dir))
    plot_generator = PlotGenerator(Path(output_dir) / "plots")

    damage_data = data_collector.collect_damage_data()
    return plot_generator.create_damage_distribution_plot(damage_data)


def create_correlation_plot(output_dir: str) -> Optional[str]:
    """
    Create a damage correlation plot.

    Args:
        output_dir: Pipeline output directory

    Returns:
        Base64 encoded plot image
    """
    data_collector = DamageDataCollector(Path(output_dir))
    plot_generator = PlotGenerator(Path(output_dir) / "plots")

    damage_data = data_collector.collect_damage_data()
    return plot_generator.create_damage_correlation_plot(damage_data)


def create_status_summary_plot(output_dir: str) -> Optional[str]:
    """
    Create a status summary plot.

    Args:
        output_dir: Pipeline output directory

    Returns:
        Base64 encoded plot image
    """
    data_collector = DamageDataCollector(Path(output_dir))
    plot_generator = PlotGenerator(Path(output_dir) / "plots")

    damage_data = data_collector.collect_damage_data()
    return plot_generator.create_status_summary_plot(damage_data)


def create_quality_damage_plot(output_dir: str) -> Optional[str]:
    """
    Create a quality vs damage plot.

    Args:
        output_dir: Pipeline output directory

    Returns:
        Base64 encoded plot image
    """
    data_collector = DamageDataCollector(Path(output_dir))
    plot_generator = PlotGenerator(Path(output_dir) / "plots")

    damage_data = data_collector.collect_damage_data()
    return plot_generator.create_quality_damage_plot(damage_data)


# Maintain backward compatibility
def get_damage_data(base_output_dir: str) -> Dict[str, Any]:
    """
    Legacy function for backward compatibility.

    Args:
        base_output_dir: Base output directory

    Returns:
        Dashboard data
    """
    return create_dashboard_data(base_output_dir)


# Export main functions for backward compatibility
__all__ = [
    "create_damage_plots",
    "create_dashboard_data",
    "create_damage_distribution_plot",
    "create_correlation_plot",
    "create_status_summary_plot",
    "create_quality_damage_plot",
    "get_damage_data",
]
