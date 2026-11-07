# Historical simulation gallery

Eight astronomical simulation recordings selected from the user-provided movie archive, reviewed 8 November 2026. Materials MD and related materials CG are excluded. Source AVI files remain untouched and are not copied into the public site.

## Regenerate

Use the project's pinned NumPy/Matplotlib stack, Pillow 11.3.0 and imageio-ffmpeg 0.6.0. From the repository root:

```sh
python scripts/build_gallery.py --source-root '/path/to/original/movie/folder'
python scripts/build_gallery.py --verify-only
python scripts/build_site.py
python scripts/check_site.py --directory dist
```

The supplied folder must contain the relative paths and source bytes recorded in selection.json. The builder checks original SHA-256 hashes, decodes complete inputs, converts to H.264/yuv420p, decodes complete outputs and compares frame counts and dimensions. PNG posters are extracted from the originals at recorded times. All selected sources were silent. The full-frame versions retain original cadence, full frame and chronology. Encoded bytes can vary with ffmpeg/platform versions; the manifest identifies the committed derivatives.

The separate 18-second edited introduction contains four 4-second excerpts and 1-second opening/closing cards. It removes the old top title frames, adds English headings and source credits, and uses hard cuts. Source imagery is fitted without upscaling, recoloring, AI restoration or generated detail. Repeated frames encode the original 3/5-fps motion in a 15-fps stream; no new simulation states or physical continuity across different clips are implied. The exact source time ranges and crops are in manifest.json.

## Quality and provenance

These are usable native-size gallery sources, rather than HD masters: 500–800 pixels wide and 3–10 fps. Granularity, limited fine detail and stepped motion are visible. Some 2D recordings have clipped lower-axis captions or incomplete visible scales; their quality notes retain those limitations. No missing labels or colors are reconstructed.

Embedded author/year/process captions are retained in full-frame versions. The caption “Koyama & Inutsuka 2002” is evidence about a recording, not a verified link to a particular paper figure or run setup. A filename's numerical suffix is not taken as an audited computational resolution. No physical or numerical validation of the simulations is claimed.

CI verifies committed media/configuration/recipe hashes and basic recorded invariants without requiring the local source archive or regenerating the movies. This is a technical media check, separate from the archive's scientific validation reports.

Two native-pixel stills are also extracted for illustration: shock-detail.png removes the old 3D title area while retaining scene axes; phase-detail.png crops a qualitative 2D field region, excluding its axes. The exact pixel rectangles and requested poster seek positions are in the manifest. Seek positions are quantized by the source cadence. No upscaling or color restoration is applied. The homepage and Ginzburg–Landau aside use these as archival illustrations, not numerical comparisons.
