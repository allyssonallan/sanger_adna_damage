"""
Report coordination for Sanger Pipeline QC reports.

This module coordinates the overall report generation process,
orchestrating statistics collection, HTML generation, and file output.
"""

import logging
from datetime import datetime
from pathlib import Path

from .statistics_collector import StatisticsCollector
from .html_template_generator import HTMLTemplateGenerator

logger = logging.getLogger(__name__)


class ReportCoordinator:
    """Coordinates the generation of comprehensive QC reports."""

    def __init__(self, output_dir: Path, damage_plotter=None):
        """
        Initialize report coordinator.

        Args:
            output_dir: Pipeline output directory
            damage_plotter: Optional damage plot generator instance
        """
        self.output_dir = Path(output_dir)
        self.report_dir = self.output_dir / "reports"
        self.report_dir.mkdir(exist_ok=True)

        # Initialize components
        self.statistics_collector = StatisticsCollector(output_dir, damage_plotter)
        self.html_generator = HTMLTemplateGenerator(output_dir)

    def generate_report(self) -> Path:
        """
        Generate comprehensive QC report.

        Returns:
            Path to generated HTML report
        """
        logger.info("Generating comprehensive QC report...")

        try:
            # Collect all statistics
            stats = self.statistics_collector.collect_pipeline_statistics()

            # Generate HTML report
            html_content = self.html_generator.generate_html_report(stats)

            # Save report
            report_file = (
                self.report_dir
                / f"sanger_qc_report_{datetime.now().strftime('%Y%m%d_%H%M%S')}.html"
            )

            with open(report_file, "w", encoding="utf-8") as f:
                f.write(html_content)

            logger.info(f"QC report generated: {report_file}")

            return report_file

        except Exception as e:
            logger.error(f"Error generating QC report: {e}")
            raise
