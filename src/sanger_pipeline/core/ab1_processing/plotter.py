"""Quality plot generation service for AB1 processing."""

from __future__ import annotations

import logging

import matplotlib

# Use a non-interactive backend so PNG generation works in headless or
# minimally configured Windows environments without Tcl/Tk.
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from Bio.SeqRecord import SeqRecord

from .models import QualityPlotRequest

logger = logging.getLogger(__name__)


class QualityPlotter:
    """Generates quality score plots from SeqRecord objects."""

    def __init__(self, min_quality: int) -> None:
        self.min_quality = min_quality

    def generate(self, record: SeqRecord, request: QualityPlotRequest) -> None:
        """Generate and save a quality plot if Phred scores are present."""
        if "phred_quality" not in record.letter_annotations:
            logger.warning("Record does not contain phred_quality annotations")
            return

        qualities = record.letter_annotations["phred_quality"]
        request.output_file.parent.mkdir(parents=True, exist_ok=True)

        plt.figure(figsize=request.figure_size)
        plt.plot(qualities, marker=".", linestyle="-", markersize=1)
        plt.axhline(
            self.min_quality,
            color="red",
            linestyle="--",
            label=f"Quality threshold ({self.min_quality})",
        )
        plt.title(f"Quality scores for {record.id}")
        plt.xlabel("Base position")
        plt.ylabel("Phred Quality")
        plt.legend()
        plt.grid(True, alpha=0.3)
        plt.tight_layout()
        plt.savefig(request.output_file, dpi=300, bbox_inches="tight")
        plt.close()
        logger.info(f"Generated quality plot: {request.output_file}")
