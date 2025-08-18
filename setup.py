"""
Setup configuration for the Sanger DNA damage analysis pipeline.
"""

from setuptools import setup, find_packages
from pathlib import Path

# Read the contents of README file
this_directory = Path(__file__).parent
long_description = (this_directory / "README.md").read_text()

setup(
    name="sanger-adna-damage",
    version="1.0.0",
    author="Allysson Allan",
    author_email="",
    description="A comprehensive pipeline for Sanger sequencing AB1 file analysis with HVS region processing and ancient DNA damage detection",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/allyssonallan/sanger_adna_damage",
    package_dir={"": "src"},
    packages=find_packages(where="src"),
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
    ],
    python_requires=">=3.8",
    install_requires=[
        "click>=8.0",
        "pyyaml>=6.0",
        "numpy>=1.20",
        "matplotlib>=3.5",
        "pandas>=1.3",
        "biopython>=1.78",  # Compatible with Python 3.8+
        "scipy>=1.7",
    ],
    extras_require={
        "dev": [
            "pytest>=7.0",
            "pytest-cov>=4.0",
            "black>=22.0",
            "flake8>=5.0",
            "mypy>=1.0",
        ],
        "docs": [
            "sphinx>=5.0",
            "sphinx-rtd-theme>=1.0",
        ],
    },
    entry_points={
        "console_scripts": [
            "sanger-pipeline=sanger_pipeline.cli.main:cli",
        ],
    },
    include_package_data=True,
    zip_safe=False,
)
