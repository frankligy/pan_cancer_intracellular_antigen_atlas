# ImmunoVerse Portal — Upgraded Edition (April 2026)

A pan-cancer atlas of therapeutic T cell targets accompanying the manuscript
*"A pan-cancer atlas of therapeutic T cell targets"* (Yarmarkovich et al., 2026).

## How to run locally
The portal is a static site. It can be opened two ways:

**Option A — Local HTTP server (recommended, works everywhere):**

```bash
cd "immuno-verse-portal"
python3 -m http.server 8000
```

Then open http://127.0.0.1:8000 in any modern browser.

**Option B — Double-click `index.html`:** works in most browsers, but Chrome
blocks `fetch()` on `file://` URLs by default. If the table stays empty, use
Option A.

## What's in this folder
- `index.html` — single-file portal (hero, atlas grid, explorer, class browser,
  manuscript highlights, methods, downloads, citation).
- `data/*.json` — 21 per-cancer tables (~38 MB total). Compact column-array
  format, pre-sorted by relative abundance.
- `data/_summary.json` — class totals, cancer metadata, labels.
- `preprocess.py` — script that regenerates the JSON files from the original
  per-cancer `*_final_enhanced_all.csv` files under
  `working files/ImmunoVerse/table/`.

## Deployment
The folder is fully self-contained and can be dropped on any static host
(Netlify, Vercel, GitHub Pages, S3, or the NYU server currently hosting
www.immuno-verse.com). No build step required.

## Data
340,916 tumor-specific HLA-presented peptides across 21 cancer types and 9
molecular aberration classes, empirically detected by Tesorai-powered
re-search of 1,823 immunopeptidomes with multi-modal safety screening.
