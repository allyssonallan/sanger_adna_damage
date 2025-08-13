"""
Helper functions for the Sanger pipeline.
"""

import logging
import yaml
from pathlib import Path
from typing import Dict, Any, Optional
from .constants import DEFAULT_CONFIG_FILE


def setup_logging(level: str = "INFO", log_file: Optional[Path] = None) -> None:
    """
    Setup logging configuration for the pipeline.

    Args:
        level: Logging level (DEBUG, INFO, WARNING, ERROR, CRITICAL)
        log_file: Optional path to log file
    """
    handlers = [logging.StreamHandler()]

    if log_file:
        handlers.append(logging.FileHandler(log_file))

    logging.basicConfig(
        level=getattr(logging, level.upper()),
        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
        handlers=handlers,
        force=True,
    )


def load_config(config_file: Optional[Path] = None) -> Dict[str, Any]:
    """
    Load configuration from YAML file.

    Args:
        config_file: Path to configuration file. If None, uses default.

    Returns:
        Configuration dictionary
    """
    if config_file is None:
        config_file = DEFAULT_CONFIG_FILE

    if not config_file.exists():
        raise FileNotFoundError(f"Configuration file not found: {config_file}")

    with open(config_file, "r") as f:
        config = yaml.safe_load(f)

    return config


def create_directories(output_dir: Path) -> Dict[str, Path]:
    """
    Create output directory structure for the pipeline.

    Args:
        output_dir: Base output directory

    Returns:
        Dictionary mapping directory names to Path objects
    """
    directories = {
        "output": output_dir,
        "fasta": output_dir / "fasta",
        "filtered": output_dir / "filtered",
        "aligned": output_dir / "aligned", 
        "consensus": output_dir / "consensus",
        "final": output_dir / "final",
        "plots": output_dir / "plots",
        "reports": output_dir / "reports"
    }

    # Create all directories
    for dir_path in directories.values():
        dir_path.mkdir(parents=True, exist_ok=True)

    return directories


def ensure_directories(directories: Dict[str, str], base_path: Path) -> Dict[str, Path]:
    """
    Create output directories if they don't exist.

    Args:
        directories: Dictionary mapping directory names to relative paths
        base_path: Base path for creating directories

    Returns:
        Dictionary mapping directory names to absolute Path objects
    """
    created_dirs = {}

    for name, rel_path in directories.items():
        dir_path = base_path / rel_path
        dir_path.mkdir(parents=True, exist_ok=True)
        created_dirs[name] = dir_path

    return created_dirs


def validate_file_exists(file_path: Path, file_type: str = "file") -> None:
    """
    Validate that a file exists and raise informative error if not.

    Args:
        file_path: Path to validate
        file_type: Type of file for error message

    Raises:
        FileNotFoundError: If file doesn't exist
    """
    if not file_path.exists():
        raise FileNotFoundError(f"{file_type.title()} not found: {file_path}")


def get_sample_name(file_path: Path, pattern: str) -> str:
    """
    Extract sample name from file path based on pattern.

    Args:
        file_path: Path to the file
        pattern: Pattern to remove from filename

    Returns:
        Sample name
    """
    name = file_path.stem
    if pattern.startswith("*"):
        suffix = pattern[1:]
        if name.endswith(suffix.replace("_filtered", "").replace(".fasta", "")):
            # Remove the pattern suffix
            if "-F" in suffix:
                name = name.replace("-F_filtered", "").replace("-F", "")
            elif "-R" in suffix:
                name = name.replace("-R_filtered", "").replace("-R", "")
            elif "HSV1" in suffix:
                name = name.replace("-HSV1_consensus", "")
            elif "HSV2" in suffix:
                name = name.replace("-HSV2_consensus", "")

    return name


def clean_dna_sequence(sequence: str) -> str:
    """
    Clean DNA sequence by removing non-ACGT characters and converting to uppercase.

    Args:
        sequence: Raw DNA sequence

    Returns:
        Cleaned DNA sequence
    """
    import re

    return re.sub(r"[^ACGTacgt]", "", sequence.upper())
