# 🔧 GitHub Pages Setup Required

## The Issue - SOLVED! ✅

The `403 permission denied` error has been fixed by updating the GitHub Actions workflow to use the modern GitHub Pages deployment method.

## ⚠️ **Action Required: Configure GitHub Pages Source**

You need to update your GitHub repository settings **once** to use the new deployment method:

### Steps:

1. **Go to your repository on GitHub**
   - Navigate to: `https://github.com/allyssonallan/sanger_adna_damage`

2. **Open Settings**
   - Click the **Settings** tab (top right of repository)

3. **Configure Pages**
   - In the left sidebar, click **Pages**
   - Under **Source**, select **GitHub Actions** (NOT "Deploy from a branch")
   - Save the settings

4. **Configure Custom Domain**
   - In the same Pages settings page
   - Under **Custom domain**, enter: `allysson.dev.br`
   - Click **Save**

## ✅ What's Fixed

### **Updated GitHub Actions Workflow**
- ✅ **Proper permissions** for Pages deployment
- ✅ **Modern actions** (`actions/deploy-pages@v2`)
- ✅ **CNAME file creation** for custom domain
- ✅ **Concurrency control** to prevent deployment conflicts

### **New Workflow Features**
- ✅ **Better error handling** and permissions
- ✅ **Automatic artifact upload** and deployment
- ✅ **Support for custom domains** with CNAME
- ✅ **Concurrent deployment protection**

## 🚀 After Configuration

Once you complete the repository settings:

1. **The workflow will run automatically** on your next push
2. **Documentation will deploy** to `https://allysson.dev.br/sanger_adna_damage/`
3. **No more permission errors** - the bot now has proper permissions

## 🔍 Monitoring Deployment

- **Actions tab**: Monitor workflow progress
- **Settings → Pages**: See deployment status and domain configuration
- **Live site**: Check `https://allysson.dev.br/sanger_adna_damage/` after ~5-10 minutes

## 📝 Technical Details

### What Changed:
```yaml
# OLD (peaceiris/actions-gh-pages)
- uses: peaceiris/actions-gh-pages@v3
  with:
    github_token: ${{ secrets.GITHUB_TOKEN }}
    publish_dir: ./docs/_build/html

# NEW (GitHub's official Pages actions)
permissions:
  pages: write
  id-token: write
  
- uses: actions/upload-pages-artifact@v2
- uses: actions/deploy-pages@v2
```

The new method uses GitHub's official Pages deployment actions which have better permission handling and reliability.
