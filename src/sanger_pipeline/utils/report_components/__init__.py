"""
Report generation components for Sanger Pipeline QC reports.

This package contains modular components for generating comprehensive
HTML reports with statistics collection, HTML template generation,
and report coordination.
"""

from .statistics_collector import StatisticsCollector
from .html_template_generator import HTMLTemplateGenerator
from .report_coordinator import ReportCoordinator

__all__ = [
    'StatisticsCollector',
    'HTMLTemplateGenerator', 
    'ReportCoordinator'
]
