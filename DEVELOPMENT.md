# Development Guide - Best Practices Applied

This document outlines the best practices implemented in this project and how to use them effectively.

## 🚀 Quick Start for Developers

```bash
# 1. Clone and enter directory
cd sanger_adna_damage

# 2. Activate virtual environment
source venv/bin/activate

# 3. Install development dependencies
make install-dev

# 4. Run tests
make test-fast

# 5. Check code quality
make quality
```

## 🧪 Testing Strategy

### Test Organization

Tests are organized following Python best practices:

``` bash
tests/
├── conftest.py                 # Shared fixtures and configuration
├── test_integration_smoke.py   # Integration/smoke tests
├── test_core/                  # Unit tests for core modules
├── test_*.py                   # Other test modules
└── originals/                  # Legacy tests (to be refactored)
```

### Test Categories

Tests are marked with pytest markers:

- `@pytest.mark.unit` - Fast unit tests (< 1 second)
- `@pytest.mark.integration` - Integration tests (1-10 seconds)
- `@pytest.mark.slow` - Slow tests (> 10 seconds)
- `@pytest.mark.requires_external` - Tests requiring external dependencies
- `@pytest.mark.requires_data` - Tests requiring test data files

### Running Tests

```bash
# Run all tests
make test

# Run only fast tests (no coverage)
make test-fast

# Run tests with coverage report
make test-cov

# Run only unit tests
pytest -m "unit"

# Run everything except slow tests
pytest -m "not slow"

# Run smoke tests only
make test-smoke

# Quick import check
make test-imports
```

## 📊 Code Quality Tools

### Linting and Formatting
```bash
# Format code with Black
make format

# Check formatting
make format-check

# Run linting
make lint

# Type checking
make type-check

# Run all quality checks
make quality
```

### Pre-commit Hooks

```bash
# Install pre-commit hooks
make pre-commit-install

# Run pre-commit on all files
make pre-commit-run
```

## 🔧 Configuration

### Project Configuration Files

1. **`pyproject.toml`** - Modern Python project configuration
   - Build system configuration
   - Dependency management
   - Tool configurations (pytest, black, mypy, coverage)

2. **`setup.cfg`** - Additional tool configuration
   - Flake8 linting rules

3. **`.pre-commit-config.yaml`** - Pre-commit hooks configuration

4. **`Makefile`** - Development task automation

### Pytest Configuration

Key settings in `pyproject.toml`:

```toml
[tool.pytest.ini_options]
testpaths = ["tests"]
markers = [
    "slow: marks tests as slow",
    "integration: marks tests as integration tests",
    "unit: marks tests as unit tests",
]
addopts = [
    "--strict-markers",
    "--cov=src/sanger_pipeline",
    "--cov-report=term-missing",
    "--cov-fail-under=80",
]
```

## 🏗️ Project Structure Best Practices

### Source Code Organization

``` bash
src/sanger_pipeline/          # Source code in src layout
├── __init__.py
├── cli/                      # Command-line interface
├── core/                     # Core business logic
├── utils/                    # Utilities and helpers
└── scripts/                  # Standalone scripts
```

### Test Organization Best Practices

1. **Mirror source structure** - Test files should mirror the source code structure
2. **Use descriptive names** - Test functions should clearly describe what they test
3. **Group related tests** - Use test classes to group related functionality
4. **Use fixtures** - Share common test setup via pytest fixtures
5. **Mark tests appropriately** - Use markers to categorize tests

## 🔄 Development Workflow

### Daily Development

1. Activate virtual environment: `source venv/bin/activate`
2. Run fast tests frequently: `make test-fast`
3. Format code before committing: `make format`
4. Run quality checks: `make quality`

### Before Committing

1. Run full test suite: `make test`
2. Check coverage: `make test-cov`
3. Run quality checks: `make quality`
4. Let pre-commit hooks run automatically

### Continuous Integration

The project includes GitHub Actions CI that:

- Tests on multiple Python versions (3.8-3.12)
- Tests on multiple OS (Ubuntu, macOS, Windows)
- Runs linting and type checking
- Generates coverage reports
- Runs smoke tests

## 📈 Coverage and Quality Metrics

### Coverage Targets

- **Unit tests**: Aim for 90%+ coverage
- **Integration tests**: Focus on critical paths
- **Overall project**: Target 80%+ coverage

### Quality Metrics

- **Linting**: All code should pass flake8
- **Type checking**: All code should pass mypy
- **Formatting**: All code should be Black-formatted

## 🛠️ Available Make Targets

Run `make help` to see all available targets:

```bash
make help                 # Show available targets
make install             # Install package in development mode
make install-dev         # Install with development dependencies
make test               # Run all tests
make test-fast          # Run tests without coverage (faster)
make test-cov           # Run tests with coverage report
make test-smoke         # Run smoke tests only
make test-imports       # Quick import check
make lint               # Run linting
make format             # Format code with black
make format-check       # Check if code is formatted correctly
make type-check         # Run type checking
make quality            # Run all quality checks
make ci                 # Run CI-like workflow
make clean              # Clean up build artifacts
make docs               # Build documentation
make build              # Build package
```

## 🌐 Browser Application Development

### In-Browser Pipeline Features

The project includes a modern web-based interface for the Sanger aDNA damage analysis pipeline. The browser application provides a complete alternative to the command-line interface.

#### Architecture

```text
browser/
├── index.html           # Main application page
├── css/
│   └── main.css        # Styling and responsive design
└── js/
    ├── main.js         # Application coordination and initialization
    ├── core.js         # Core utilities and data management
    ├── pipeline.js     # Processing pipeline implementation
    ├── ui-manager.js   # User interface and interactions
    ├── abif-parser.js  # ABIF file format parser
    ├── chromatogram.js # Chromatogram visualization
    ├── sequence-processor.js # Sequence processing algorithms
    └── export-manager.js     # Data export functionality
```

#### File Loading System (🚧 In Development)

The browser application features multiple file loading methods to handle various use cases and browser limitations:

##### 1. **📁 Traditional File Picker**

- Standard HTML5 file input with multiple selection
- Accept filter for `.ab1` files
- Comprehensive error handling and validation
- Compatible with all modern browsers

##### 2. **📂 Drag & Drop Interface**

- Large visual drop zone with hover feedback
- Automatic file type filtering
- Visual animations and user feedback
- Click-to-open file picker functionality

##### 3. **🗂️ File System Access API**

- Modern browser API (Chrome/Edge 86+)
- Direct file system access with proper permissions
- Enhanced file selection with native OS dialogs
- Fallback detection for unsupported browsers

##### 4. **📂 HTTP Server Loading**

- Direct loading from project's `/input` directory
- Modal selection interface for choosing specific files
- Local development server with CORS support
- Bypasses browser file security restrictions

##### 5. **🎯 Demo Data Generator**

- Synthetic sequence data for testing
- Bypasses file system entirely
- Immediate testing capability
- Pre-configured sample data sets

##### 6. **📋 Paste Support**

- Automatic clipboard file detection
- Background processing without UI interaction
- Supports files copied from file managers
- Cross-platform compatibility

#### Testing and Debugging Tools

##### File System Diagnostics

```javascript
// Test file loading system
app.testFileLoading();

// Load demo data for testing
app.loadDemoData();

// Open File System API dialog
app.openFileSystemDialog();
```

##### Development Server

```bash
# Start local development server
python3 scripts/serve_files.py

# Access browser application
open http://localhost:8080/browser/
```

#### Browser Compatibility

| Feature | Chrome | Firefox | Safari | Edge |
|---------|---------|---------|---------|-------|
| File Picker | ✅ | ✅ | ✅ | ✅ |
| Drag & Drop | ✅ | ✅ | ✅ | ✅ |
| File System API | ✅ | ❌ | ❌ | ✅ |
| Paste Support | ✅ | ✅ | ⚠️ | ✅ |
| HTTP Loading | ✅ | ✅ | ✅ | ✅ |

#### Current Development Status

- ✅ **Core Pipeline**: Complete ES6 modular architecture
- ✅ **File Parsing**: ABIF format parser with trace data extraction
- ✅ **UI Components**: Responsive design with theme support
- ✅ **Data Export**: ZIP archive generation with multiple formats
- 🚧 **File Loading**: Multiple loading methods (in testing)
- 🚧 **Chromatogram Display**: Canvas-based visualization
- 📋 **Advanced Features**: Damage analysis visualization (planned)

#### Known Issues and Limitations

1. **File System Access**: Limited to Chromium-based browsers
2. **Large Files**: Memory constraints for very large AB1 files
3. **CORS Restrictions**: Local file access requires development server
4. **Mobile Support**: Touch interactions need optimization

#### Future Enhancements

- **WebAssembly Integration**: For performance-critical operations
- **Progressive Web App**: Offline capability and app-like experience
- **Real-time Collaboration**: Multi-user analysis sessions
- **Cloud Storage Integration**: Direct loading from cloud services

## 🔍 Troubleshooting

### Common Issues

1. **Import Errors**: Ensure virtual environment is activated and dependencies installed
2. **Test Failures**: Check if external dependencies are available
3. **Coverage Issues**: Coverage below 80% will fail CI - focus on testing critical paths
4. **Linting Errors**: Run `make format` to auto-fix many issues

### Test-Specific Issues

1. **Missing fixtures**: Check `conftest.py` for available fixtures
2. **Slow tests**: Use `pytest -m "not slow"` to skip slow tests during development
3. **External dependencies**: Mark tests with `@pytest.mark.requires_external`

## 📚 Best Practices Summary

### ✅ What We Implemented

1. **Proper test organization** - Tests in dedicated directory with logical structure
2. **Pytest configuration** - Markers, coverage, and proper test discovery
3. **Code quality tools** - Black, flake8, mypy, pre-commit hooks
4. **Development automation** - Makefile for common tasks
5. **CI/CD setup** - GitHub Actions for automated testing
6. **Modern Python packaging** - pyproject.toml configuration
7. **Test categorization** - Markers for different test types
8. **Coverage reporting** - HTML and terminal coverage reports

### ❌ What We Fixed

1. **Tests outside test directory** - Moved functionality to proper test structure
2. **CLI command mismatch** - Fixed test expectations to match actual CLI
3. **Missing test configuration** - Added comprehensive pytest configuration
4. **No code quality enforcement** - Added linting and formatting tools
5. **Manual test running** - Automated with make targets and scripts

### 🎯 Key Benefits

1. **Faster development** - Quick feedback from fast tests
2. **Higher code quality** - Automated linting and formatting
3. **Better maintainability** - Clear test organization and documentation
4. **Easier collaboration** - Consistent development practices
5. **Continuous integration** - Automated testing on multiple platforms

This setup provides a solid foundation for professional Python development with excellent testing practices and code quality enforcement.
