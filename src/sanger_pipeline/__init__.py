"""
Sanger DNA Damage Analysis Pipeline

A comprehensive pipeline for processing Sanger sequencing AB1 files,
including quality control, alignment, consensus building, and damage analysis.
"""

__version__ = "1.0.0"
__author__ = "Allysson Allan"

from .core.pipeline import SangerPipeline
from .core.ab1_converter import AB1Converter
from .core.consensus_builder import ConsensusBuilder
from .core.quality_filter import QualityFilter
from .core.adna_damage_analyzer import ADNADamageAnalyzer

__all__ = [
    "SangerPipeline",
    "AB1Converter",
    "ConsensusBuilder",
    "QualityFilter",
    "ADNADamageAnalyzer",
]
