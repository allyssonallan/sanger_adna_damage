=============
Contributing
=============

Thank you for your interest in contributing to the Sanger DNA Damage Analysis Pipeline! This guide will help you get started with contributing code, documentation, or other improvements.

🤝 Ways to Contribute
=====================

There are many ways to contribute to this project:

* **Bug Reports**: Help us identify and fix issues
* **Feature Requests**: Suggest new functionality or improvements
* **Code Contributions**: Submit bug fixes, new features, or optimizations
* **Documentation**: Improve existing docs or add new guides
* **Testing**: Help test the pipeline with different data types
* **Community Support**: Help other users in discussions and issues

📋 Getting Started
==================

Prerequisites
-------------

Before contributing, ensure you have:

* Python 3.8+ installed
* Git version control system
* Basic familiarity with the pipeline (complete the :doc:`quickstart` guide)
* Development tools (IDE/editor of your choice)

Setting Up Development Environment
---------------------------------

1. **Fork the Repository**

   Click the "Fork" button on the GitHub repository page to create your own copy.

2. **Clone Your Fork**

   .. code-block:: bash

      git clone https://github.com/YOUR_USERNAME/sanger_adna_damage.git
      cd sanger_adna_damage

3. **Set Up Remote**

   .. code-block:: bash

      # Add upstream remote
      git remote add upstream https://github.com/allyssonallan/sanger_adna_damage.git
      
      # Verify remotes
      git remote -v

4. **Create Development Environment**

   .. code-block:: bash

      # Create virtual environment
      python3 -m venv dev_env
      source dev_env/bin/activate  # On Windows: dev_env\\Scripts\\activate
      
      # Install development dependencies
      pip install -r requirements.txt
      pip install -e .
      
      # Install development tools (if available)
      pip install -r requirements-dev.txt

5. **Verify Installation**

   .. code-block:: bash

      # Test the pipeline
      python -m src.sanger_pipeline.cli.main --help
      
      # Run tests (if available)
      pytest tests/

🐛 Bug Reports
==============

Found a Bug?
------------

Before creating a bug report:

1. **Check existing issues** to avoid duplicates
2. **Try the latest version** to see if it's already fixed
3. **Follow troubleshooting guide** to rule out common issues

Creating a Bug Report
---------------------

Use this template for bug reports:

.. code-block:: markdown

   ## Bug Report
   
   **Description**
   A clear and concise description of what the bug is.
   
   **To Reproduce**
   Steps to reproduce the behavior:
   1. Run command: `python -m src.sanger_pipeline.cli.main ...`
   2. With input files: `...`
   3. See error: `...`
   
   **Expected Behavior**
   What you expected to happen.
   
   **Screenshots/Logs**
   If applicable, add error messages or log outputs.
   
   **Environment**
   - OS: [e.g. macOS 12.0, Ubuntu 20.04]
   - Python version: [e.g. 3.9.7]
   - Pipeline version: [e.g. 1.0.0]
   - MAFFT version: [e.g. 7.490]
   
   **Additional Context**
   Any other context about the problem.

💡 Feature Requests
===================

Suggesting New Features
----------------------

Feature requests are welcome! Before submitting:

1. **Check existing requests** to avoid duplicates
2. **Consider the scope** - does it fit the pipeline's goals?
3. **Think about implementation** - is it technically feasible?

Feature Request Template
-----------------------

.. code-block:: markdown

   ## Feature Request
   
   **Is your feature request related to a problem?**
   A clear description of what the problem is.
   
   **Describe the solution you'd like**
   A clear description of what you want to happen.
   
   **Describe alternatives you've considered**
   Alternative solutions or features you've considered.
   
   **Use Cases**
   Specific examples of how this feature would be used.
   
   **Additional Context**
   Any other context, mockups, or examples.

💻 Code Contributions
====================

Development Workflow
--------------------

1. **Create a Branch**

   .. code-block:: bash

      # Sync with upstream
      git fetch upstream
      git checkout main
      git merge upstream/main
      
      # Create feature branch
      git checkout -b feature/my-new-feature

2. **Make Changes**

   * Follow the coding standards (see below)
   * Write tests for new functionality
   * Update documentation as needed
   * Commit changes with clear messages

3. **Test Your Changes**

   .. code-block:: bash

      # Run tests
      pytest tests/
      
      # Test with sample data
      python -m src.sanger_pipeline.cli.main run-pipeline \
          --input-dir ./test_data \
          --output-dir ./test_output

4. **Push and Create Pull Request**

   .. code-block:: bash

      # Push to your fork
      git push origin feature/my-new-feature
      
      # Create pull request on GitHub

Coding Standards
---------------

**Python Style**:
* Follow PEP 8 style guide
* Use meaningful variable and function names
* Add docstrings to all public functions and classes
* Keep functions focused and concise

**Example Function**:

.. code-block:: python

   def calculate_damage_score(sequences: List[str], positions: int = 20) -> float:
       """Calculate ancient DNA damage score from sequences.
       
       Args:
           sequences: List of DNA sequences to analyze
           positions: Number of positions to analyze from each end
           
       Returns:
           Damage score between 0 and 1
           
       Raises:
           ValueError: If sequences list is empty or positions < 1
       """
       if not sequences:
           raise ValueError("Sequences list cannot be empty")
       
       if positions < 1:
           raise ValueError("Positions must be >= 1")
       
       # Implementation here
       return damage_score

**Documentation Style**:
* Use Google-style docstrings
* Include type hints for function parameters and returns
* Document all parameters, return values, and exceptions

**Testing Standards**:
* Write unit tests for all new functions
* Include integration tests for new features
* Test error conditions and edge cases
* Aim for >80% code coverage

**Example Test**:

.. code-block:: python

   import pytest
   from src.sanger_pipeline.utils.damage_analyzer import calculate_damage_score

   def test_calculate_damage_score_valid_input():
       """Test damage score calculation with valid sequences."""
       sequences = ["ATCGATCG", "TTCGATCA", "ATCGATCG"]
       score = calculate_damage_score(sequences)
       assert 0.0 <= score <= 1.0

   def test_calculate_damage_score_empty_input():
       """Test that empty sequences raise ValueError."""
       with pytest.raises(ValueError, match="Sequences list cannot be empty"):
           calculate_damage_score([])

Commit Message Guidelines
------------------------

Write clear, descriptive commit messages:

.. code-block:: text

   Short (50 chars or less) summary of changes

   More detailed explanatory text, if necessary. Wrap it to about 72
   characters. The blank line separating the summary from the body is
   critical (unless you omit the body entirely).

   Further paragraphs come after blank lines.

   - Bullet points are okay, too
   - Use a hyphen or asterisk for the bullet

**Examples**:

.. code-block:: text

   Add bootstrap analysis for damage assessment
   
   Fix AB1 conversion error with corrupted files
   
   Update documentation for new HVS region feature
   
   Improve performance of consensus building algorithm

📝 Documentation Contributions
==============================

Types of Documentation
----------------------

* **User Guides**: Help users accomplish specific tasks
* **Tutorials**: Step-by-step learning experiences
* **API Documentation**: Technical reference for developers
* **How-To Guides**: Solutions to common problems

Documentation Standards
-----------------------

* Use clear, concise language
* Include code examples for technical content
* Test all code examples to ensure they work
* Use reStructuredText (.rst) format for Sphinx
* Follow the established documentation structure

Building Documentation Locally
------------------------------

.. code-block:: bash

   # Install documentation dependencies
   pip install sphinx sphinx_rtd_theme
   
   # Build documentation
   cd docs/
   make html
   
   # View documentation
   open _build/html/index.html  # macOS
   xdg-open _build/html/index.html  # Linux

🧪 Testing Guidelines
====================

Test Categories
---------------

**Unit Tests**:
* Test individual functions in isolation
* Mock external dependencies
* Fast execution (<1 second per test)

**Integration Tests**:
* Test component interactions
* Use real dependencies where appropriate
* Medium execution time (1-10 seconds per test)

**End-to-End Tests**:
* Test complete workflows
* Use real AB1 files (small test dataset)
* Slower execution (10+ seconds per test)

Running Tests
------------

.. code-block:: bash

   # Run all tests
   pytest
   
   # Run specific test file
   pytest tests/test_damage_analyzer.py
   
   # Run with coverage
   pytest --cov=src
   
   # Run only fast tests
   pytest -m "not slow"

Test Data
--------

* Use small, synthetic test files when possible
* Include a few real AB1 files for integration testing
* Document the source and characteristics of test data
* Keep test data in the `tests/data/` directory

🔄 Pull Request Process
=======================

Pull Request Checklist
----------------------

Before submitting a pull request:

- [ ] Code follows the project's coding standards
- [ ] All tests pass locally
- [ ] New code has appropriate test coverage
- [ ] Documentation is updated (if applicable)
- [ ] Commit messages are clear and descriptive
- [ ] Branch is up-to-date with main branch

Pull Request Template
---------------------

.. code-block:: markdown

   ## Description
   Brief description of what this PR does.
   
   ## Type of Change
   - [ ] Bug fix (non-breaking change which fixes an issue)
   - [ ] New feature (non-breaking change which adds functionality)
   - [ ] Breaking change (fix or feature that would cause existing functionality to not work as expected)
   - [ ] Documentation update
   
   ## Testing
   - [ ] I have added tests that prove my fix is effective or that my feature works
   - [ ] New and existing unit tests pass locally with my changes
   - [ ] I have tested this with real AB1 files
   
   ## Checklist
   - [ ] My code follows the style guidelines of this project
   - [ ] I have performed a self-review of my own code
   - [ ] I have commented my code, particularly in hard-to-understand areas
   - [ ] I have made corresponding changes to the documentation
   - [ ] My changes generate no new warnings
   
   ## Screenshots/Examples
   If applicable, add examples of the changes.

Review Process
-------------

1. **Automated Checks**: CI/CD will run tests and checks
2. **Code Review**: Maintainers will review your code
3. **Feedback**: Address any requested changes
4. **Approval**: Once approved, your PR will be merged

🏷️ Release Process
==================

Version Numbering
-----------------

We follow Semantic Versioning (SemVer):

* **MAJOR**: Breaking changes (e.g., 1.0.0 → 2.0.0)
* **MINOR**: New features, backward compatible (e.g., 1.0.0 → 1.1.0)
* **PATCH**: Bug fixes, backward compatible (e.g., 1.0.0 → 1.0.1)

Changelog
--------

We maintain a changelog following the "Keep a Changelog" format:

.. code-block:: markdown

   # Changelog
   
   ## [Unreleased]
   ### Added
   - New bootstrap analysis feature
   
   ### Changed
   - Improved performance of consensus building
   
   ### Fixed
   - Fixed AB1 conversion bug with certain file types
   
   ## [1.0.0] - 2024-01-15
   ### Added
   - Initial release of the pipeline
   - Complete AB1 to consensus workflow
   - Ancient DNA damage analysis

🌟 Recognition
==============

Contributors
-----------

We recognize all contributors in:

* README.md contributors section
* Documentation acknowledgments
* Release notes
* Git commit history

Types of Recognition
-------------------

* **Code Contributors**: Direct code contributions
* **Bug Reporters**: High-quality bug reports
* **Documentation Contributors**: Documentation improvements
* **Community Contributors**: Helping others, discussions

📞 Communication
===============

Where to Discuss
----------------

* **GitHub Issues**: Bug reports, feature requests
* **GitHub Discussions**: General questions, ideas, help
* **Pull Requests**: Code review and discussion

Communication Guidelines
------------------------

* Be respectful and inclusive
* Use clear, concise language
* Provide context and examples
* Search existing discussions before posting
* Stay on topic

Code of Conduct
---------------

We are committed to providing a welcoming and inclusive environment. All contributors are expected to:

* Use welcoming and inclusive language
* Be respectful of differing viewpoints and experiences
* Gracefully accept constructive criticism
* Focus on what is best for the community
* Show empathy towards other community members

🎯 Getting Started with Your First Contribution
===============================================

Good First Issues
----------------

Look for issues labeled:
* `good first issue`: Easy problems to get started
* `help wanted`: Issues where we'd appreciate help
* `documentation`: Documentation improvements needed

Simple Contribution Ideas
-------------------------

1. **Fix Typos**: Documentation or code comments
2. **Add Examples**: More usage examples in documentation
3. **Improve Error Messages**: Make error messages more helpful
4. **Add Tests**: Increase test coverage
5. **Performance Improvements**: Optimize slow operations

Steps for First Contribution
----------------------------

1. **Find an Issue**: Look through open issues for something interesting
2. **Comment**: Let maintainers know you're working on it
3. **Ask Questions**: Don't hesitate to ask for clarification
4. **Start Small**: Begin with a small, focused change
5. **Learn from Feedback**: Use code review as a learning opportunity

Thank you for contributing to the Sanger DNA Damage Analysis Pipeline! Your contributions help make this tool better for the entire research community.
