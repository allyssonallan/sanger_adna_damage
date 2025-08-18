"""
Test configuration for the Sanger pipeline tests.
"""

import sys
import pytest
from pathlib import Path
import tempfile
import shutil
import os

# Make repository root importable so `import src.sanger_pipeline` works
REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT))


@pytest.fixture(scope="session")
def repo_root():
    """Repository root directory."""
    return REPO_ROOT


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
        "damage": {
            "threshold": 0.05,
            "bootstrap_iterations": 1000,
            "significance_level": 0.05,
        },
        "bootstrap": {"iterations": 1000, "confidence_level": 0.95},
        "output": {
            "directories": {
                "fasta": "fasta",
                "filtered": "filtered",
                "consensus": "consensus",
                "plots": "plots",
                "aligned": "aligned",
                "final": "final",
                "logs": "logs",
            }
        },
    }


@pytest.fixture
def test_input_dir(temp_dir):
    """Create a test input directory structure."""
    input_dir = temp_dir / "input"
    input_dir.mkdir(parents=True, exist_ok=True)
    return input_dir


@pytest.fixture
def test_output_dir(temp_dir):
    """Create a test output directory structure."""
    output_dir = temp_dir / "output"
    output_dir.mkdir(parents=True, exist_ok=True)
    return output_dir


@pytest.fixture
def sample_fasta_content():
    """Sample FASTA content for testing."""
    return """
>sample1_HVS1_F
ATCGATCGATCGATCG
>sample1_HVS1_R
CGATCGATCGATCGAT
>sample2_HVS2_F
GCTAGCTAGCTAGCTA
>sample2_HVS2_R
TAGCTAGCTAGCTAG
""".strip()


@pytest.fixture
def mock_config_file(temp_dir, sample_config):
    """Create a mock configuration file."""
    import yaml

    config_file = temp_dir / "test_config.yaml"
    with open(config_file, "w") as f:
        yaml.dump(sample_config, f)

    return config_file


@pytest.fixture
def clean_environment():
    """Ensure clean environment for tests."""
    # Store original environment
    original_env = os.environ.copy()

    # Clean up any pipeline-specific environment variables
    env_vars_to_clean = ["SANGER_CONFIG", "SANGER_DEBUG", "SANGER_LOG_LEVEL"]

    for var in env_vars_to_clean:
        if var in os.environ:
            del os.environ[var]

    yield

    # Restore original environment
    os.environ.clear()
    os.environ.update(original_env)


# Pytest configuration for test selection
def pytest_configure(config):
    """Configure pytest with custom markers."""
    config.addinivalue_line("markers", "unit: marks tests as unit tests")
    config.addinivalue_line("markers", "integration: marks tests as integration tests")
    config.addinivalue_line(
        "markers", "slow: marks tests as slow (deselect with '-m \"not slow\"')"
    )
    config.addinivalue_line(
        "markers", "requires_external: marks tests that require external dependencies"
    )
    config.addinivalue_line(
        "markers", "requires_data: marks tests that require test data files"
    )


def pytest_collection_modifyitems(config, items):
    """Modify test collection to add markers automatically."""
    for item in items:
        # Add 'unit' marker to tests that don't have integration/slow markers
        if not any(
            marker.name in ["integration", "slow", "requires_external"]
            for marker in item.iter_markers()
        ):
            item.add_marker(pytest.mark.unit)

        # Add 'slow' marker to tests with 'integration' marker
        if item.get_closest_marker("integration"):
            item.add_marker(pytest.mark.slow)
