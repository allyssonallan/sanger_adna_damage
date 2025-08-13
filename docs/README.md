# Documentation for Sanger DNA Damage Analysis Pipeline

This directory contains the complete Sphinx documentation for the Sanger DNA Damage Analysis Pipeline.

## Quick Start

### Prerequisites

```bash
# Install documentation dependencies
pip install sphinx sphinx_rtd_theme sphinx-autobuild myst-parser
```

### Building Documentation

```bash
# Build HTML documentation
make html

# View documentation
open _build/html/index.html  # macOS
xdg-open _build/html/index.html  # Linux
```

### Development

```bash
# Install development dependencies
make dev-install

# Build with error checking
make dev-build

# Auto-rebuild on changes (development server)
make livehtml
```

## Documentation Structure

``` bash
docs/
├── source/                    # Source files
│   ├── index.rst             # Main documentation index
│   ├── installation.rst      # Installation guide
│   ├── quickstart.rst        # Quick start guide
│   ├── configuration.rst     # Configuration reference
│   ├── cli_reference.rst     # CLI command reference
│   ├── understanding_damage_analysis.rst  # Damage analysis theory
│   ├── troubleshooting.rst   # Troubleshooting guide
│   ├── contributing.rst      # Contributing guidelines
│   ├── changelog.rst         # Version history
│   ├── license.rst           # License information
│   ├── tutorials/            # Tutorial directory
│   │   ├── index.rst         # Tutorial index
│   │   └── first_analysis.rst # First analysis tutorial
│   ├── howto/                # How-to guides
│   │   ├── index.rst         # How-to index
│   │   └── process_single_sample.rst # Single sample guide
│   ├── api/                  # API documentation
│   │   ├── index.rst         # API index
│   │   └── core.rst          # Core module documentation
│   ├── _static/              # Static files (CSS, JS, images)
│   ├── _templates/           # Custom templates
│   └── conf.py               # Sphinx configuration
├── Makefile                  # Build commands
└── README.md                 # This file
```

## Documentation Types

### User Documentation

- **Installation Guide**: Step-by-step installation instructions
- **Quick Start**: Get up and running in 5 minutes
- **Configuration**: Complete configuration reference
- **CLI Reference**: All command-line options and examples
- **Tutorials**: Step-by-step learning guides
- **How-To Guides**: Solution-focused guides for specific tasks

### Scientific Documentation

- **Understanding Damage Analysis**: Theory and interpretation
- **Ancient DNA Workflows**: Specialized guides for aDNA research
- **Statistical Methods**: Bootstrap analysis and significance testing

### Developer Documentation

- **API Reference**: Complete module and class documentation
- **Contributing Guide**: How to contribute to the project
- **Architecture**: System design and component interaction

### Support Documentation

- **Troubleshooting**: Common issues and solutions
- **FAQ**: Frequently asked questions
- **Changelog**: Version history and updates

## Building Different Formats

```bash
# HTML (default)
make html

# PDF (requires LaTeX)
make latexpdf

# Check links
make linkcheck

# Clean build directory
make clean

# Full rebuild
make rebuild
```

## Live Development

For documentation development with auto-rebuild:

```bash
# Start development server
make livehtml

# Navigate to http://localhost:8000
# Documentation rebuilds automatically on file changes
```

## Writing Documentation

### Style Guide

- Use clear, concise language
- Include code examples for all procedures
- Test all code examples to ensure they work
- Use consistent terminology throughout
- Follow the established structure and formatting

### reStructuredText Basics

```rst
Chapter Title
=============

Section Title
-------------

Subsection Title
~~~~~~~~~~~~~~~~

**Bold text**
*Italic text*
``Code text``

.. code-block:: bash

   # Code block example
   python -m src.sanger_pipeline.cli.main --help

.. note::
   This is a note box.

.. warning::
   This is a warning box.
```

### Cross-References

```rst
# Link to other documents
See :doc:`installation` for setup instructions.

# Link to sections
See :ref:`configuration-basics` for basic setup.

# Link to API
See :class:`src.sanger_pipeline.core.pipeline.SangerPipeline`
```

## API Documentation

API documentation is automatically generated from docstrings:

```bash
# Generate API docs from source code
make apidoc

# Build with API docs included
make rebuild
```

### Docstring Format

Use Google-style docstrings:

```python
def process_sample(input_file: str, output_dir: str) -> bool:
    """Process a single AB1 sample file.
    
    Args:
        input_file: Path to the AB1 input file
        output_dir: Directory for output files
        
    Returns:
        True if processing succeeded, False otherwise
        
    Raises:
        FileNotFoundError: If input file doesn't exist
        PermissionError: If output directory isn't writable
        
    Example:
        >>> success = process_sample("sample.ab1", "./output")
        >>> print(f"Processing {'succeeded' if success else 'failed'}")
    """
```

## Quality Checks

```bash
# Check documentation for errors
make check

# Check external links
make linkcheck

# Check spelling (if enabled)
make spelling
```

## Contributing to Documentation

1. **Find areas to improve**: Look for unclear sections or missing content
2. **Follow the style guide**: Maintain consistency with existing documentation
3. **Test your examples**: Ensure all code examples work correctly
4. **Build locally**: Test your changes before submitting
5. **Submit pull request**: Include a clear description of your changes

### Common Improvements Needed

- More step-by-step tutorials
- Additional how-to guides for specific use cases
- Better code examples and real-world scenarios
- Improved troubleshooting coverage
- More detailed API documentation

## Documentation Themes

The documentation uses the Read the Docs theme with custom styling:

- **Responsive design**: Works on desktop and mobile
- **Dark/light mode**: User-selectable theme
- **Search functionality**: Full-text search across all documents
- **Navigation**: Sidebar navigation with collapsible sections
- **Interactive elements**: Tabs, accordions, and code blocks

## Deployment

Documentation is automatically built and deployed when:

- Changes are pushed to the main branch
- Pull requests are merged
- Releases are tagged

Local builds can be deployed manually:

```bash
# Build for production
make html

# Deploy to GitHub Pages (if configured)
# Copy _build/html/* to gh-pages branch
```

## Getting Help

If you need help with documentation:

1. **Check existing docs**: Look for similar examples
2. **Sphinx documentation**: https://www.sphinx-doc.org/
3. **reStructuredText guide**: https://docutils.sourceforge.io/rst.html
4. **GitHub issues**: Report documentation bugs or requests
5. **Community discussions**: Ask questions in GitHub Discussions

## Maintenance

Regular documentation maintenance tasks:

- **Update screenshots**: Keep visual elements current
- **Check links**: Ensure external links still work
- **Review accuracy**: Verify technical content is up-to-date
- **Add new features**: Document new functionality as it's added
- **Improve clarity**: Continuously improve explanations based on user feedback

The documentation is a living resource that grows and improves with the project and community feedback.
