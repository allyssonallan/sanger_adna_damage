============
Installation
============

This guide walks you through installing the Sanger DNA Damage Analysis Pipeline on your system.

.. note::
   **About This Tool**
   
   Before installing, please note that this pipeline is designed for **sample prioritization and preliminary screening**, not for definitive ancient DNA authentication. The tool helps identify promising samples for follow-up NGS analysis based on damage patterns, insert size, and haplogroup assessment.

📋 Prerequisites
================

System Requirements
-------------------

* **Python**: 3.8 or higher
* **Operating System**: Linux, macOS, or Windows
* **Memory**: Minimum 4GB RAM (8GB+ recommended for large datasets)
* **Storage**: 1GB free space for installation + space for your data

Required External Tools
-----------------------

The pipeline requires several external bioinformatics tools:

MAFFT (Multiple Sequence Alignment)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Ubuntu/Debian:**

.. code-block:: bash

   sudo apt update
   sudo apt install mafft

**macOS (using Homebrew):**

.. code-block:: bash

   brew install mafft

**Windows:**

Download from `MAFFT website <https://mafft.cbrc.jp/alignment/software/>`_ and follow installation instructions.

**Verify Installation:**

.. code-block:: bash

   mafft --version

Bio Module Dependencies
~~~~~~~~~~~~~~~~~~~~~~~

The pipeline uses BioPython which should be installed automatically, but you may need additional system packages:

**Ubuntu/Debian:**

.. code-block:: bash

   sudo apt install python3-dev python3-pip

**macOS:**

.. code-block:: bash

   # Usually pre-installed with Python

🚀 Installation Methods
=======================

Method 1: Development Installation (Recommended)
------------------------------------------------

This method is recommended for users who want to modify or contribute to the pipeline.

1. **Clone the Repository**

   .. code-block:: bash

      git clone https://github.com/allyssonallan/sanger_adna_damage.git
      cd sanger_adna_damage

2. **Create Virtual Environment**

   .. code-block:: bash

      # Create virtual environment
      python3 -m venv venv
      
      # Activate virtual environment
      # On Linux/macOS:
      source venv/bin/activate
      
      # On Windows:
      venv\\Scripts\\activate

3. **Install Dependencies**

   .. code-block:: bash

      # Upgrade pip
      pip install --upgrade pip
      
      # Install Python dependencies
      pip install -r requirements.txt
      
      # Install the package in development mode
      pip install -e .

4. **Verify Installation**

   .. code-block:: bash

      # Test the CLI
      python -m src.sanger_pipeline.cli.main --help
      
      # Test import
      python -c "from src.sanger_pipeline.core.pipeline import SangerPipeline; print('Installation successful!')"

Method 2: Direct Installation from Source
-----------------------------------------

For users who want a cleaner installation without development files:

1. **Download and Extract**

   .. code-block:: bash

      # Download the latest release
      wget https://github.com/allyssonallan/sanger_adna_damage/archive/main.zip
      unzip main.zip
      cd sanger_adna_damage-main

2. **Install**

   .. code-block:: bash

      pip install -r requirements.txt
      pip install .

🔧 Configuration Setup
======================

1. **Copy Default Configuration**

   .. code-block:: bash

      # Copy the default configuration to a working directory
      cp config/default_config.yaml my_config.yaml

2. **Edit Configuration** (Optional)

   .. code-block:: bash

      # Edit with your preferred editor
      nano my_config.yaml

3. **Verify Configuration**

   .. code-block:: bash

      python -m src.sanger_pipeline.cli.main status --config my_config.yaml

📁 Directory Structure Setup
============================

Create your working directories:

.. code-block:: bash

   # Create project structure
   mkdir -p my_sanger_project/{input,output,config}
   
   # Copy configuration
   cp config/default_config.yaml my_sanger_project/config/
   
   # Your directory structure should look like:
   my_sanger_project/
   ├── input/          # Place your AB1 files here
   ├── output/         # Pipeline results will go here
   └── config/         # Configuration files

🧪 Testing Installation
=======================

Test with Sample Data
---------------------

.. code-block:: bash

   # Create test directories
   mkdir -p test_run/{input,output}
   
   # If you have sample AB1 files, place them in test_run/input/
   # Then run a test analysis:
   
   python -m src.sanger_pipeline.cli.main run-pipeline \\
       --input-dir ./test_run/input \\
       --output-dir ./test_run/output \\
       --config ./config/default_config.yaml

Run Unit Tests
--------------

.. code-block:: bash

   # Run the test suite (if available)
   python -m pytest tests/

Verify All Components
---------------------

.. code-block:: bash

   # Check external tools
   mafft --version
   
   # Check Python modules
   python -c "import Bio; print(f'BioPython version: {Bio.__version__}')"
   
   # Check pipeline modules
   python -c "from src.sanger_pipeline.core.pipeline import SangerPipeline; print('Pipeline import successful')"
   
   # Generate help to verify CLI
   python -m src.sanger_pipeline.cli.main --help

🔍 Troubleshooting Installation
===============================

Common Issues
-------------

**ImportError: No module named 'Bio'**

.. code-block:: bash

   pip install biopython

**MAFFT command not found**

Make sure MAFFT is installed and in your PATH:

.. code-block:: bash

   which mafft  # Should show the path to mafft
   echo $PATH   # Check if mafft directory is in PATH

**Permission denied errors**

Use virtual environments or install with --user flag:

.. code-block:: bash

   pip install --user -r requirements.txt

**Python version issues**

Check your Python version:

.. code-block:: bash

   python --version  # Should be 3.8+

Virtual Environment Issues
--------------------------

If you have problems with virtual environments:

.. code-block:: bash

   # Remove existing environment
   rm -rf venv
   
   # Create new environment with explicit Python version
   python3.8 -m venv venv  # or python3.9, python3.10, etc.
   
   source venv/bin/activate
   pip install --upgrade pip
   pip install -r requirements.txt

💡 Development Setup
====================

For developers who want to contribute:

1. **Fork the Repository** on GitHub

2. **Clone Your Fork**

   .. code-block:: bash

      git clone https://github.com/YOUR_USERNAME/sanger_adna_damage.git
      cd sanger_adna_damage

3. **Set Up Development Environment**

   .. code-block:: bash

      # Create development environment
      python3 -m venv dev_env
      source dev_env/bin/activate
      
      # Install development dependencies
      pip install -r requirements.txt
      pip install -e .
      
      # Install development tools (if requirements-dev.txt exists)
      pip install -r requirements-dev.txt

4. **Set Up Pre-commit Hooks** (if available)

   .. code-block:: bash

      pre-commit install

🎯 Next Steps
=============

Once installation is complete:

1. **Read the Quick Start Guide**: :doc:`quickstart`
2. **Configure the Pipeline**: :doc:`configuration`
3. **Follow Tutorials**: :doc:`tutorials/index`
4. **Run Your First Analysis**: :doc:`howto/index`

If you encounter any issues during installation, please check the :doc:`troubleshooting` guide or open an issue on GitHub.
