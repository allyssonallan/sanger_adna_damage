.PHONY: help install install-dev test test-fast test-cov lint format clean docs build
.DEFAULT_GOAL := help

# Variables
PYTHON := python
VENV := venv
SRC := src
TESTS := tests

help: ## Show this help message
	@echo "Available targets:"
	@grep -E '^[a-zA-Z_-]+:.*?## .*$$' $(MAKEFILE_LIST) | sort | awk 'BEGIN {FS = ":.*?## "}; {printf "  \033[36m%-15s\033[0m %s\n", $$1, $$2}'

install: ## Install package in development mode
	$(PYTHON) -m pip install -e .

install-dev: ## Install package with development dependencies
	$(PYTHON) -m pip install -e ".[dev,test]"

test: ## Run all tests
	$(PYTHON) -m pytest $(TESTS)/ -v

test-fast: ## Run tests without coverage (faster)
	$(PYTHON) -m pytest $(TESTS)/ -v --no-cov

test-cov: ## Run tests with coverage report
	$(PYTHON) -m pytest $(TESTS)/ -v --cov=$(SRC) --cov-report=html --cov-report=term-missing

test-smoke: ## Run smoke tests only
	$(PYTHON) run_tests.py --mode smoke

test-imports: ## Quick import check
	$(PYTHON) run_tests.py --mode imports

lint: ## Run linting (flake8)
	$(PYTHON) -m flake8 $(SRC) $(TESTS)

lint-report: ## Generate detailed linting report
	$(PYTHON) -m flake8 $(SRC) $(TESTS) --count --statistics --tee --output-file=flake8-report.txt

lint-critical: ## Check only critical linting errors
	$(PYTHON) -m flake8 $(SRC) $(TESTS) --select=E9,F63,F7,F82 --show-source

lint-fix: ## Auto-fix common linting issues
	@echo "Fixing trailing whitespace and blank lines..."
	$(PYTHON) -c "import os, re; [open(f, 'w').write(re.sub(r'[ \t]+$$', '', open(f).read(), flags=re.MULTILINE)) for root, dirs, files in os.walk('src') for f in [os.path.join(root, file) for file in files if file.endswith('.py')] if os.path.isfile(f)]"
	$(PYTHON) -c "import os, re; [open(f, 'w').write(re.sub(r'[ \t]+$$', '', open(f).read(), flags=re.MULTILINE)) for root, dirs, files in os.walk('tests') for f in [os.path.join(root, file) for file in files if file.endswith('.py')] if os.path.isfile(f)]"
	@echo "Running black to fix formatting..."
	$(PYTHON) -m black $(SRC) $(TESTS) run_tests.py setup.py

format: ## Format code with black
	$(PYTHON) -m black $(SRC) $(TESTS) run_tests.py setup.py

format-check: ## Check if code is formatted correctly
	$(PYTHON) -m black --check $(SRC) $(TESTS) run_tests.py setup.py

type-check: ## Run type checking with mypy
	$(PYTHON) -m mypy $(SRC)

clean: ## Clean up build artifacts and cache
	rm -rf build/
	rm -rf dist/
	rm -rf *.egg-info/
	rm -rf htmlcov/
	rm -rf .coverage
	rm -rf .pytest_cache/
	rm -rf .mypy_cache/
	find . -type d -name __pycache__ -exec rm -rf {} +
	find . -type f -name "*.pyc" -delete

docs: ## Build documentation
	cd docs && make html

docs-serve: ## Serve documentation locally
	cd docs/_build/html && $(PYTHON) -m http.server 8000

build: ## Build package
	$(PYTHON) -m build

# Quality checks
quality: lint type-check format-check ## Run all quality checks

# CI-like workflow
ci: install-dev quality test-cov ## Run CI-like workflow

# Development setup
setup-dev: ## Set up development environment
	$(PYTHON) -m venv $(VENV)
	$(VENV)/bin/pip install --upgrade pip setuptools wheel
	$(VENV)/bin/pip install -e ".[dev,test]"
	@echo "Development environment set up in $(VENV)/"
	@echo "Activate with: source $(VENV)/bin/activate"

# Pre-commit hooks setup
pre-commit-install: ## Install pre-commit hooks
	$(PYTHON) -m pre_commit install

pre-commit-run: ## Run pre-commit hooks on all files
	$(PYTHON) -m pre_commit run --all-files
