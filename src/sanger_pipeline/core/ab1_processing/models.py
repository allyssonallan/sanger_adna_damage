"""Data models used by AB1 processing services."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Tuple


@dataclass(frozen=True)
class AB1ConversionRequest:
    """Input parameters for converting an AB1 file to FASTA."""

    ab1_file: Path
    output_file: Path


@dataclass(frozen=True)
class QualityPlotRequest:
    """Input parameters for generating a quality plot."""

    output_file: Path
    figure_size: Tuple[int, int] = (12, 4)
