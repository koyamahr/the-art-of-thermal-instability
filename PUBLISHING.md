# Publication rehearsal and automatic updates

## Current status

Initial publication succeeded on 2026-11-03. Site: https://koyamahr.github.io/the-art-of-thermal-instability/ . Repository: https://github.com/koyamahr/the-art-of-thermal-instability . `PAGES_AUTO_DEPLOY=true` is enabled: pushes to main validate and deploy automatically.

The first [rehearsal](https://github.com/koyamahr/the-art-of-thermal-instability/actions/runs/41558943576) and [manual deployment](https://github.com/koyamahr/the-art-of-thermal-instability/actions/runs/41559093976) succeeded, including all 49 scientific checks. The steps below document how to repeat the rehearsal and initial setup.

## 1. Rehearse locally

Use the project environment described in `experiments/field-1965/README.md`, then run:

```sh
.venv/bin/python experiments/field-1965/publish.py
.venv/bin/python experiments/fgh-1969/publish.py
python3 scripts/check_site.py
python3 scripts/build_site.py
python3 scripts/check_site.py --directory dist
python3 -m http.server 8000 --bind 127.0.0.1 --directory dist
```

Visit http://127.0.0.1:8000 and inspect all three exhibits and their downloads. Stop the server with Ctrl+C. Rebuilding Field and FGH/GHF currently runs 60 scientific checks. This step does not publish to the internet.

## 2. Connect GitHub and rehearse there

The repository `koyamahr/the-art-of-thermal-instability` has been created. Authenticate Git on this Mac, then push this project with `git push -u origin main`. For a free personal GitHub plan, use a public repository for Pages. Author PDFs are ignored by Git; only the checked `dist/` artifact, generated from public `docs/` inputs, is packaged for the website.

For an initial setup on another repository, do not add the repository variable `PAGES_AUTO_DEPLOY` yet. On this repository it is already enabled. With that variable absent, pushing `main` runs checks and creates the `site-rehearsal` artifact without deploying. In Actions, open **Check, rehearse, and publish site**, inspect the logs, and download the artifact to review it with a local HTTP server.

For another rehearsal, choose **Run workflow**, select **main**, leave **publish** unchecked, and run. A successful rehearsal verifies the remote calculation/build process, but does not produce a public preview URL.

## 3. First real publication

In repository **Settings → Pages → Build and deployment → Source**, select **GitHub Actions**. This replaces the earlier branch-and-`/docs` publishing instructions.

In **Actions → Check, rehearse, and publish site → Run workflow**, choose **main** and check **publish**. A successful run deploys only the generated `dist/` artifact to Pages. Open the URL returned by the deployment and check the home page, Figures 1–3, and code/CSV download links. Resolve any Pages environment approval or account-specific repository setting shown by GitHub.

The normal project URL is `https://<account>.github.io/the-art-of-thermal-instability/`; use the actual deployment URL as authoritative. Publishing from a non-main branch is explicitly rejected. Original implementation code is MIT; content/data/media have separate terms. See LICENSING.md.

## 4. Enable automatic updates

This repository already has the variable below. For a new setup, after the first successful publication, add **Settings → Secrets and variables → Actions → Variables → New repository variable**:

- Name: `PAGES_AUTO_DEPLOY`
- Value: `true`

Then each push to `main` recomputes the figures, checks the site, and deploys it if every step passes. Pull requests run checks and save a rehearsal artifact, but never deploy. The workflow uses GitHub's built-in token; no personal access token is required by the workflow.

Automatic updates mean publishing changes committed and pushed to GitHub. Saving a local file alone does not update the public site. There is no scheduled content generation or automatic literature collection.

## Fixes and rollback

If validation fails, deployment does not start and the existing site remains in place. Fix the change and push again. To restore earlier content, revert the unwanted commit on `main` and push; when automatic publishing is enabled, the reverted version is rebuilt and deployed. Removing the variable disables future automatic deployments without deleting the existing site.

## Runtime and files

The GitHub workflow uses Ubuntu and Python 3.11 with pinned direct dependencies. This differs from the initial local Python 3.9.6 run; each regenerated validation report records its actual runtime. `environment.txt` is the original local dependency snapshot, not a claim about the CI environment. The workflow does not push generated files back into Git: they are part of the deployment artifact, while experiment sources remain the source of truth.

Official references: [custom Pages workflows](https://docs.github.com/en/pages/getting-started-with-github-pages/using-custom-workflows-with-github-pages), [publishing source configuration](https://docs.github.com/en/pages/getting-started-with-github-pages/configuring-a-publishing-source-for-your-github-pages-site).

## Searchable archive build

After scientific outputs and reading paths are prepared, `scripts/build_site.py` creates all 47 paper overviews at `papers/<id>/`, embeds the complete catalog in the homepage, and adds page-specific canonical and Open Graph metadata to all 65 HTML pages. Existing exhibit paths and catalog IDs stay unchanged. The generated sitemap lists every HTML page. `scripts/check_site.py --directory dist` verifies the static catalog, paper scopes, search metadata, sitemap, local links and scientific provenance before either upload step. The build uses only Python standard-library modules and never reads private source PDFs.
