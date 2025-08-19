"""
CLI command modules for the Sanger pipeline.

This package contains modular command groups organized by functionality.
"""

from .pipeline_commands import run_pipeline, status, pipeline
from .conversion_commands import (
    convert_ab1,
    convert_to_hsd,
    convert_ab1_enhanced,
    validate_primers,
    generate_primer_config,
)
from .analysis_commands import analyze_damage, damage_summary
from .hsd_commands import hsd
from .report_commands import generate_report

__all__ = [
    "run_pipeline",
    "status",
    "convert_ab1",
    "convert_ab1_enhanced",
    "convert_to_hsd",
    "validate_primers",
    "generate_primer_config",
    "analyze_damage",
    "damage_summary",
    "hsd",
    "generate_report",
]
