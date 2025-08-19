"""
Sanger DNA Damage Analysis Pipeline

A pipeline for analyzing ancient DNA damage patterns in Sanger sequencing data.
"""

__version__ = "0.1.0"

# Core classes for external use
from .core.enhanced_ab1_converter_fixed import EnhancedAB1Converter as AB1Converter
from .core.pipeline import SangerPipeline

__all__ = [
    "AB1Converter",  # Now points to EnhancedAB1Converter for backward compatibility
    "SangerPipeline",
]
