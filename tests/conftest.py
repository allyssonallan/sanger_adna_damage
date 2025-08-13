"""
Test configuration for the Sanger pipeline tests.
"""

import pytest
from pathlib import Path
import tempfile
import shutil

@pytest.fixture
def temp_dir():
    """Create a temporary directory for tests."""
    temp_dir = Path(tempfile.mkdtemp())
    yield temp_dir
    shutil.rmtree(temp_dir)

@pytest.fixture
def sample_config():
    """Sample configuration for testing."""
    return {
        "quality": {"min_phred_score": 20},
        "alignment": {"tool": "mafft", "parameters": "--auto"},
        "output": {
            "directories": {
                "fasta": "fasta",
                "filtered": "filtered",
                "consensus": "consensus",
                "plots": "plots",
                "aligned": "aligned",
                "final": "final",
                "logs": "logs"
            }
        }
    }
