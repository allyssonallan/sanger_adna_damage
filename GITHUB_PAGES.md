# GitHub Pages Deployment Guide

This repository supports multiple ways to deploy documentation to GitHub Pages.

## Option 1: Automatic CI/CD Deployment (Recommended)

The repository includes a GitHub Actions workflow that automatically builds and deploys documentation when you push to main.

### Setup:
1. Go to your repository settings on GitHub
2. Navigate to **Settings** → **Pages**
3. Under **Source**, select **GitHub Actions**
4. The workflow will automatically deploy on the next push to `main`

### Workflow file: `.github/workflows/docs.yml`
- Builds documentation using Sphinx
- Deploys to `gh-pages` branch automatically
- Runs on every push to `main`

## Option 2: Manual GitHub Pages from `/gh-pages` folder

### Build and deploy locally:
```bash
cd docs
./build_docs.sh
git add ../gh-pages/
git commit -m "Update documentation"
git push
```

### Setup GitHub Pages:
1. Go to **Settings** → **Pages**
2. Under **Source**, select **Deploy from a branch**
3. Choose **main** branch and **/gh-pages** folder
4. Save settings

## Option 3: Deploy to separate `gh-pages` branch

### One-time setup:
```bash
# Create and switch to gh-pages branch
git checkout --orphan gh-pages
git rm -rf .
echo "GitHub Pages branch" > README.md
git add README.md
git commit -m "Initial gh-pages commit"
git push origin gh-pages
git checkout main
```

### Deploy docs:
```bash
cd docs
make html
cd ..
git checkout gh-pages
cp -r docs/_build/html/* .
git add .
git commit -m "Update documentation"
git push origin gh-pages
git checkout main
```

## Current Setup

This repository is configured for **Option 1** (CI/CD) and **Option 2** (gh-pages folder) simultaneously:

- ✅ GitHub Actions workflow ready in `.github/workflows/docs.yml`
- ✅ Build script creates `gh-pages/` directory
- ✅ `.gitignore` allows committing `gh-pages/` content
- ✅ README.md links updated to GitHub Pages URLs

## Documentation URLs

Once deployed, your documentation will be available at:
- **Main documentation**: https://allyssonallan.github.io/sanger_adna_damage/
- **Installation guide**: https://allyssonallan.github.io/sanger_adna_damage/installation.html
- **Quick start**: https://allyssonallan.github.io/sanger_adna_damage/quickstart.html
- **API reference**: https://allyssonallan.github.io/sanger_adna_damage/api/

## Troubleshooting

### Documentation not updating
1. Check GitHub Actions tab for build errors
2. Verify Pages settings point to correct source
3. Check that `gh-pages/` directory contains HTML files

### 404 errors
1. Ensure `index.html` exists in the deployment location
2. Check that GitHub Pages is enabled
3. Verify the repository is public (or you have GitHub Pro for private repos)

### Build failures
1. Check all dependencies are listed in `requirements.txt`
2. Verify Sphinx configuration in `docs/source/conf.py`
3. Run `cd docs && ./build_docs.sh` locally to test
