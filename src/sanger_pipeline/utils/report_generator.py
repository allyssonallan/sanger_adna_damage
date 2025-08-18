"""
Comprehensive QC Report Generator for Sanger Pipeline.

This module provides a simplified interface to the modular report generation
system, maintaining backward compatibility while using the new component-based
architecture.
"""

import logging
from pathlib import Path

from .damage_plots import DamagePlotGenerator
from .report_components import ReportCoordinator

logger = logging.getLogger(__name__)


class QCReportGenerator:
    """Generate comprehensive QC reports for Sanger pipeline analysis."""
    
    def __init__(self, output_dir: Path):
        """
        Initialize report generator.
        
        Args:
            output_dir: Pipeline output directory
        """
        self.output_dir = Path(output_dir)
        
        # Initialize damage plot generator
        self.damage_plotter = DamagePlotGenerator(output_dir)
        
        # Initialize report coordinator with damage plotter
        self.report_coordinator = ReportCoordinator(output_dir, self.damage_plotter)
    
    def collect_pipeline_statistics(self):
        """
        Collect pipeline statistics.
        
        Returns:
            Statistics dictionary
        """
        return self.report_coordinator.statistics_collector.collect_pipeline_statistics()
    
    def generate_html_report(self, stats):
        """
        Generate HTML report from statistics.
        
        Args:
            stats: Statistics dictionary
            
        Returns:
            HTML content string
        """
        return self.report_coordinator.html_generator.generate_html_report(stats)
    
    def generate_report(self) -> Path:
        """
        Generate comprehensive QC report.
        
        Returns:
            Path to generated HTML report
        """
        return self.report_coordinator.generate_report()
