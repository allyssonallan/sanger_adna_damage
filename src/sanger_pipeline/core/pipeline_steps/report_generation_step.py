"""
Report generation step handler for the Sanger pipeline.

This module handles comprehensive QC report generation.
"""

import logging
from pathlib import Path
from typing import Dict

logger = logging.getLogger(__name__)


class ReportGenerationStep:
    """Handles comprehensive QC report generation."""
    
    def execute(self, directories: Dict[str, Path]) -> Dict:
        """
        Execute report generation step.
        
        Args:
            directories: Pipeline output directories
            
        Returns:
            Dictionary with generation status
        """
        logger.info("Step 5: Generating comprehensive QC report")

        # Include damage analysis results in the report
        damage_dir = directories["output"] / "damage_analysis"
        if damage_dir.exists():
            logger.info("Including aDNA damage analysis results in report")

        # Generate comprehensive QC report
        try:
            from ...utils.report_generator import QCReportGenerator
            
            report_generator = QCReportGenerator(directories["output"])
            report_file = report_generator.generate_report()
            
            logger.info(f"QC report generated successfully: {report_file}")
            return {"success": True, "report_path": str(report_file)}
            
        except Exception as e:
            logger.error(f"Failed to generate QC report: {e}")
            logger.info("Report generation failed, but pipeline continues")
            return {"success": False, "error": str(e)}
