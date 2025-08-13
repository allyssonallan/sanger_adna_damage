#!/bin/bash
# Documentation build and deployment script for GitHub Pages

set -e  # Exit on any error

echo "🚀 Building Sanger DNA Damage Analysis Pipeline Documentation"
echo "============================================================="

# Check if we're in the docs directory
if [ ! -f "source/conf.py" ]; then
    echo "❌ Error: Must be run from the docs/ directory"
    echo "Usage: cd docs && ./build_docs.sh"
    exit 1
fi

# Activate virtual environment if it exists
if [ -d "../venv" ]; then
    echo "🐍 Activating virtual environment..."
    source ../venv/bin/activate
    echo "✅ Virtual environment activated"
else
    echo "⚠️  No virtual environment found, using system Python"
fi

# Check for required dependencies
echo "📋 Checking dependencies..."

check_dependency() {
    if ! python -c "import $1" 2>/dev/null; then
        echo "❌ Missing dependency: $1"
        echo "Install with: pip install $1"
        exit 1
    else
        echo "✅ $1 available"
    fi
}

check_dependency "sphinx"
check_dependency "furo"

# Optional dependencies
if python -c "import myst_parser" 2>/dev/null; then
    echo "✅ myst_parser available (optional)"
else
    echo "⚠️  myst_parser not available (optional for Markdown support)"
fi

echo ""

# Clean previous builds
echo "🧹 Cleaning previous builds..."
if [ -d "_build" ]; then
    rm -rf _build
    echo "✅ Cleaned _build directory"
fi

# Build HTML documentation
echo "🔨 Building HTML documentation..."
make html

if [ $? -eq 0 ]; then
    echo "✅ HTML documentation built successfully"
else
    echo "❌ HTML build failed"
    exit 1
fi

# Check for build warnings
if [ -f "_build/html/.doctrees/environment.pickle" ]; then
    echo "✅ Build completed without critical errors"
else
    echo "⚠️  Build may have had warnings"
fi

# Validate the build
echo ""
echo "🔍 Validating documentation..."

# Check if index.html exists
if [ -f "_build/html/index.html" ]; then
    echo "✅ Main index page created"
else
    echo "❌ Main index page missing"
    exit 1
fi

# Check if key pages exist
key_pages=(
    "installation.html"
    "quickstart.html"
    "configuration.html"
    "cli_reference.html"
    "understanding_damage_analysis.html"
    "troubleshooting.html"
)

for page in "${key_pages[@]}"; do
    if [ -f "_build/html/$page" ]; then
        echo "✅ $page exists"
    else
        echo "❌ $page missing"
        exit 1
    fi
done

# Check tutorials and howto directories
if [ -d "_build/html/tutorials" ]; then
    echo "✅ Tutorials directory exists"
else
    echo "⚠️  Tutorials directory missing"
fi

if [ -d "_build/html/howto" ]; then
    echo "✅ How-to guides directory exists"
else
    echo "⚠️  How-to guides directory missing"
fi

# Check API documentation
if [ -d "_build/html/api" ]; then
    echo "✅ API documentation directory exists"
else
    echo "⚠️  API documentation directory missing"
fi

# Get build statistics
echo ""
echo "📊 Build Statistics:"
html_files=$(find _build/html -name "*.html" | wc -l)
echo "   📄 HTML files: $html_files"

css_files=$(find _build/html -name "*.css" | wc -l)
echo "   🎨 CSS files: $css_files"

js_files=$(find _build/html -name "*.js" | wc -l)
echo "   ⚡ JavaScript files: $js_files"

# Calculate total size
total_size=$(du -sh _build/html | cut -f1)
echo "   💾 Total size: $total_size"

# Optional: Run link checker
echo ""
echo "🔗 Checking links (optional)..."
if make linkcheck 2>/dev/null; then
    echo "✅ Link check completed (see _build/linkcheck/output.txt for results)"
else
    echo "⚠️  Link check skipped (may require internet connection)"
fi

# Success message
echo ""
echo "🎉 Documentation build completed successfully!"
echo ""
echo "📖 View documentation:"
echo "   Local file: file://$(pwd)/_build/html/index.html"
echo ""
echo "🚀 To serve locally:"
echo "   cd _build/html && python -m http.server 8000"
echo "   Then open: http://localhost:8000"
echo ""

# Deploy to GitHub Pages directory
echo "📄 Deploying to GitHub Pages..."
if [ -d "../gh-pages" ]; then
    echo "🗑️  Cleaning existing gh-pages directory..."
    rm -rf ../gh-pages/*
else
    echo "📁 Creating gh-pages directory..."
    mkdir -p ../gh-pages
fi

echo "📋 Copying built documentation to gh-pages..."
cp -r _build/html/* ../gh-pages/

echo "✅ Documentation deployed to gh-pages/ directory"
echo "💡 To publish on GitHub:"
echo "   1. Commit and push the gh-pages/ directory"
echo "   2. Go to GitHub → Settings → Pages"
echo "   3. Select 'Deploy from a branch'"
echo "   4. Choose 'main' branch and '/gh-pages' folder"
echo ""

# Optional: Open in browser (macOS/Linux)
if command -v open >/dev/null 2>&1; then
    read -p "🌐 Open documentation in browser? [y/N]: " -n 1 -r
    echo
    if [[ $REPLY =~ ^[Yy]$ ]]; then
        open "_build/html/index.html"
    fi
elif command -v xdg-open >/dev/null 2>&1; then
    read -p "🌐 Open documentation in browser? [y/N]: " -n 1 -r
    echo
    if [[ $REPLY =~ ^[Yy]$ ]]; then
        xdg-open "_build/html/index.html"
    fi
fi

echo "✨ Build and deployment script completed!"
