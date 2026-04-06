"""AB1 and FASTA I/O service."""

from __future__ import annotations

import logging
from pathlib import Path

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

from ...utils.helpers import validate_file_exists
from .models import AB1ConversionRequest

logger = logging.getLogger(__name__)


class AB1IO:
    """Handles AB1 read and FASTA write operations."""

    def convert_to_fasta(self, request: AB1ConversionRequest) -> SeqRecord:
        """Read AB1 input and write a FASTA output with stable naming."""
        validate_file_exists(request.ab1_file, "AB1 file")

        try:
            record = SeqIO.read(request.ab1_file, "abi")
            logger.info(f"Read AB1 file: {request.ab1_file}")
        except Exception as exc:
            raise ValueError(f"Failed to parse AB1 file {request.ab1_file}: {exc}")

        record.id = request.output_file.name
        record.description = ""
        self.write_fasta(record, request.output_file)
        logger.info(f"Wrote FASTA file: {request.output_file}")
        return record

    def write_fasta(self, record: SeqRecord, output_file: Path) -> None:
        """Write a SeqRecord as FASTA."""
        output_file.parent.mkdir(parents=True, exist_ok=True)
        SeqIO.write(record, output_file, "fasta")
