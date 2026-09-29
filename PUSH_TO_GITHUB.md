# Push to GitHub — Instructions

Your local git repo is ready at `/home/user/hydrate_app` (branch `main`, 2 commits).

## Option A: Create new repo on GitHub UI (recommended, manual)

1. Go to **https://github.com/new**
2. Repository name: `hydrate-formation-app` (or your choice)
   - Description: `Clathrate hydrate equilibria — PRSV2 + vdW-P model, interactive web app`
   - Visibility: Public
   - **Do NOT** initialize with README, .gitignore, or license (we already have them)
3. Click **Create repository**
4. GitHub will show you the remote URL, e.g.:
   ```
   https://github.com/YOUR_USERNAME/hydrate-formation-app.git
   ```

5. On your local machine (after downloading and unzipping `hydrate_app.zip`):

```bash
cd hydrate_app
git init
git config user.email "qazi@qazinasir.com"
git config user.name "Qazi Nasir"
git add .
git commit -m "Initial commit: hydrate equilibria app v8"
git branch -M main
git remote add origin https://github.com/YOUR_USERNAME/hydrate-formation-app.git
git push -u origin main
```

If your repo was already `git init`ed (this zip includes .git? No, the zip above excludes .git), use:

```bash
cd hydrate_app
git remote add origin https://github.com/YOUR_USERNAME/hydrate-formation-app.git
git branch -M main
git push -u origin main
```

## Option B: Let this environment push for you (needs PAT)

If you want me to push directly from this sandbox:

1. Create the empty repo on GitHub as in Option A step 1-3
2. Create a Personal Access Token:
   - Go to https://github.com/settings/tokens/new
   - Select **classic token**, expiration 7 days, scope `repo`
   - Generate and copy token (starts with `ghp_...`)
3. Paste here: `https://github.com/YOUR_USERNAME/hydrate-formation-app.git` and the token

Then I will run:

```bash
cd /home/user/hydrate_app
git remote add origin https://<TOKEN>@github.com/YOUR_USERNAME/hydrate-formation-app.git
git push -u origin main
```

After push, remove token from remote:

```bash
git remote set-url origin https://github.com/YOUR_USERNAME/hydrate-formation-app.git
```

## Current local state

```
branch: main
commits:
16b8795 Add copyright Qazi Nasir to README
c9e3f46 Initial commit: hydrate equilibria app v8 with Qazi Nasir copyright

Files: 77 (zip excludes .git and *_BAD.*)
```

## GitHub Pages (optional)

To host the app as a website:
- Repo Settings → Pages → Source: Deploy from branch `main`, folder `/ (root)`
- Your app will be at `https://YOUR_USERNAME.github.io/hydrate-formation-app/`
- Since `index.html` is self-contained, it works immediately on Pages.

© 2026 Qazi Nasir — https://qazinasir.com
