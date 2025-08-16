# GitHub Pages Deployment Guide

This repository supports multiple ways to deploy documentation to GitHub Pages using your custom domain `allysson.dev.br`.

## 🌐 Custom Domain Configuration

Your GitHub account is configured with a custom domain (`allysson.dev.br`), which means:
- ✅ **Your docs URL**: `https://allysson.dev.br/sanger_adna_damage/`
- ❌ **Not available**: `https://allyssonallan.github.io/sanger_adna_damage/`

## Option 1: Automatic CI/CD Deployment (Recommended)

The repository includes a GitHub Actions workflow that automatically builds and deploys documentation when you push to main.

### Setup

1. Go to your repository settings on GitHub
2. Navigate to **Settings** → **Pages**
3. Under **Source**, select **GitHub Actions**
4. The workflow will automatically deploy on the next push to `main`

### Workflow file: `.github/workflows/docs.yml`

- Builds documentation using Sphinx
- Deploys to `gh-pages` branch automatically
- Runs on every push to `main`

## Option 2: Manual GitHub Pages from `/gh-pages` folder

### Build and deploy locally

```bash
cd docs
./build_docs.sh
git add ../gh-pages/
git commit -m "Update documentation"
git push
```

### Setup GitHub Pages

1. Go to **Settings** → **Pages**
2. Under **Source**, select **Deploy from a branch**
3. Choose **main** branch and **/gh-pages** folder
4. Save settings

## 🔧 Alternative Domain Options

### Option A: Keep Using Your Custom Domain (Recommended)

**Pros:**
- ✅ Professional appearance with your domain
- ✅ Consistent with your other projects
- ✅ Better for branding and portfolio

**Cons:**
- ❌ All repos use the same domain
- ❌ Cannot use separate `github.io` subdomain

### Option B: Create GitHub Organization for Separate Domain

If you want a dedicated domain like `sanger-adna-tools.github.io`:

1. Create a new GitHub organization (e.g., `sanger-adna-tools`)
2. Transfer the repository to the organization
3. Configure GitHub Pages (will get its own `.github.io` subdomain)

### Option C: Remove Custom Domain (Not Recommended)

To use `allyssonallan.github.io` instead:

1. Go to your **Personal Settings** → **Pages**
2. Remove the custom domain configuration
3. **Warning**: This affects ALL your repositories' GitHub Pages

## Current Setup

This repository is configured for **Option 1** (CI/CD) with your custom domain:

- ✅ GitHub Actions workflow ready in `.github/workflows/docs.yml`
- ✅ Build script creates `gh-pages/` directory
- ✅ `.gitignore` allows committing `gh-pages/` content
- ✅ README.md links updated to custom domain URLs

## 📍 Documentation URLs

Your documentation will be available at:
- **Main documentation**: https://allysson.dev.br/sanger_adna_damage/
- **Installation guide**: https://allysson.dev.br/sanger_adna_damage/installation.html
- **Quick start**: https://allysson.dev.br/sanger_adna_damage/quickstart.html
- **API reference**: https://allysson.dev.br/sanger_adna_damage/api/

## 🐛 Troubleshooting

### Documentation not updating

1. Check GitHub Actions tab for build errors
2. Verify Pages settings point to correct source
3. Check that `gh-pages/` directory contains HTML files

### 404 errors

1. Ensure `index.html` exists in the deployment location
2. Check that GitHub Pages is enabled
3. Verify the repository is public (or you have GitHub Pro for private repos)
4. **Custom domain**: Ensure your domain DNS is properly configured

### Build failures

1. Check all dependencies are listed in `requirements.txt`
2. Verify Sphinx configuration in `docs/source/conf.py`
3. Run `cd docs && ./build_docs.sh` locally to test

### Custom Domain Issues

1. Check that your domain's DNS CNAME points to `allyssonallan.github.io`
2. Verify SSL certificate is active
3. Test domain resolution: `nslookup allysson.dev.br`

## 💡 Recommendation

**Keep using your custom domain** (`allysson.dev.br/sanger_adna_damage/`). It looks more professional and is consistent with your other projects. The GitHub Actions workflow is already configured to work with your custom domain setup.
