# ihc_method

A [workflowr][] project.

[workflowr]: https://github.com/workflowr/workflowr

## Deploying the site to Netlify

`docs/` holds the rendered site and is committed, so Netlify has nothing to build —
it only serves a directory. `netlify.toml` at the repo root already says so:

```toml
[build]
  publish = "docs"
  command = ""
```

`publish = "docs"` is what makes `docs/index.html` the landing page: Netlify serves
`index.html` from the publish root, so `https://<site>/` resolves to it. Pointing
Netlify at the repository root instead serves a directory with no `index.html`, and
every visitor gets a 404.

**One-time setup.** On [app.netlify.com](https://app.netlify.com): *Add new site →
Import an existing project → GitHub → `sceriff0/ihc_method`*. The build settings are
read from `netlify.toml`, so leave them as offered — publish directory `docs`, build
command empty. Netlify then redeploys on every push to `main`.

**The normal loop** is unchanged; deployment is a side effect of pushing:

```r
renv::restore()                          # first time only
workflowr::wflow_build("analysis/*.Rmd") # or wflow_publish() to tie HTML to a commit
```
```sh
git add docs && git commit -m ":rocket: rebuild site" && git push
```

Netlify's build image has no R, no renv and no pandoc — that is why `command` is
empty. Knitting happens on your machine and git carries the HTML.

**Drag-and-drop alternative**, if you want a preview without connecting the repo:
drop the `docs/` folder on [app.netlify.com/drop](https://app.netlify.com/drop). Same
result, but nothing redeploys on push.

Two things to know:

- **Knit before the first deploy.** Only `index`, `about` and `license` are currently
  in `docs/`; every analysis page needs a build (and the off-repo data in `data/`)
  before it exists to serve.
- **`docs/.nojekyll` is a GitHub Pages artefact.** Harmless on Netlify — leave it if
  you also publish via Pages, since Pages needs it to serve `site_libs/`.
