"""
Constants and configuration values for the Sanger pipeline.
"""

from pathlib import Path

# File extensions
AB1_EXTENSION = ".ab1"
FASTA_EXTENSION = ".fasta"
PNG_EXTENSION = ".png"

# Default quality thresholds
DEFAULT_MIN_QUALITY = 20
DEFAULT_MIN_SEQUENCE_LENGTH = 30
DEFAULT_TERMINAL_LENGTH = 10
DEFAULT_BOOTSTRAP_ITERATIONS = 10000

# Default directories
DEFAULT_OUTPUT_DIRS = {
    "fasta": "fasta",
    "filtered": "filtered",
    "consensus": "consensus",
    "plots": "plots",
    "aligned": "aligned",
    "final": "final",
    "logs": "logs",
}

# File patterns
PATTERNS = {
    "ab1_files": "*.ab1",
    "forward_filtered": "*-F_filtered.fasta",
    "reverse_filtered": "*-R_filtered.fasta",
    "hvs1_consensus": "*-HVS1_consensus.fasta",
    "hvs2_consensus": "*-HVS2_consensus.fasta",
    "hvs3_consensus": "*-HVS3_consensus.fasta",
}

# Configuration file locations
PROJECT_ROOT = Path(__file__).parent.parent.parent
CONFIG_DIR = PROJECT_ROOT / "config"
DEFAULT_CONFIG_FILE = CONFIG_DIR / "default_config.yaml"

# Environment variables
PYTHON_ENV_VAR = "PYTHON"
QUARTO_ENV_VAR = "QUARTO"
DEFAULT_PYTHON = "python3"
DEFAULT_QUARTO = "quarto"
