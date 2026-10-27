#!/usr/bin/env python3
"""Shared results I/O so the batch and the GUI write ONE canonical result set,
and a manual re-save of a single field updates the aggregate.

Canonical layout under <outdir> (e.g. results_all/):
  percell/<stem>.csv    per-cell rows for one field (SOURCE OF TRUTH)
  masks/<stem>.masks.tif
  overlays/<stem>.png
  all_cells.csv             = concat of percell/*.csv        (rebuilt)
  summary_per_field.csv     one row per field                (rebuilt)
  summary_per_condition.csv line x treatment means           (rebuilt)

Re-saving one field overwrites its percell/mask/overlay, then rebuild_aggregate()
regenerates all_cells + summaries — so the corrected field replaces its old rows.
"""
import glob
import os
import re

import numpy as np
import pandas as pd
import tifffile
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from skimage.segmentation import find_boundaries
from skimage.morphology import binary_dilation, disk

# metric columns summarised (only those present are used)
_METRICS = ["membrane_mean", "bright_membrane_mean", "cytoplasm_mean",
            "whole_cell_mean", "membrane_to_cyto_ratio", "membrane_frac",
            "bright_membrane_cover", "cytoplasm_frac", "puncta_frac",
            "puncta_cover_fixed", "puncta_cover_adaptive", "puncta_intensity_sum"]


def parse_meta(name):
    """(line, treatment, field) from a filename/stem."""
    b = os.path.basename(name)
    line = ("R172H" if b.startswith("R172H")
            else "TKCC10" if b.startswith("TKCC") else "syn43")
    tl = b.lower()
    treat = ("Veh" if "veh" in tl
             else "PF-846" if ("pf846" in tl or "pf-846" in tl)
             else "progesterone" if "progester" in tl else "unknown")
    mf = re.search(r"(r\d+c\d+f\d+)", b)
    return line, treat, (mf.group(1) if mf else "?")


def save_overlay(path, img, masks, res, cover_note=""):
    vmax = np.percentile(img, 99.5)
    base = np.clip(img / vmax, 0, 1)
    fig, ax = plt.subplots(1, 5, figsize=(50, 10))

    # 1) raw Filipin image
    ax[0].imshow(base, cmap="gray"); ax[0].set_title("Filipin (raw)", fontsize=15)

    # 2) segmented cells (yellow outlines)
    ov = np.dstack([base] * 3)
    ov[binary_dilation(find_boundaries(masks, mode="outer"), disk(2))] = [1, 1, 0]
    ax[1].imshow(ov); ax[1].set_title(f"{int(masks.max())} cells", fontsize=15)

    # 3) ring membrane (green) + cytoplasm (red) — the intensity compartments
    seg = np.zeros((*masks.shape, 3))
    seg[..., 0][res["cytoplasm"]] = 0.35
    seg[res["membrane"]] = [0, 0.9, 0]
    ax[2].imshow(base, cmap="gray"); ax[2].imshow(seg, alpha=0.6)
    ax[2].set_title("ring membrane (green) / cytoplasm (red)", fontsize=13)

    # 4) genuinely bright 'real' membrane (magenta), dim green ring for reference
    bm = np.zeros((*masks.shape, 3))
    bm[res["membrane"]] = [0, 0.4, 0]
    bm[res["bright_membrane"]] = [1, 0, 1]
    bcov = f"{100*res['df'].bright_membrane_cover.mean():.0f}% of ring"
    ax[3].imshow(base, cmap="gray"); ax[3].imshow(bm, alpha=0.75)
    ax[3].set_title(f"bright membrane (magenta) — {bcov}", fontsize=13)

    # 5) puncta (cyan)
    pov = np.dstack([base] * 3)
    pov[binary_dilation(res["pmask"], disk(1))] = [0, 1, 1]
    ax[4].imshow(pov); ax[4].set_title(f"puncta (fixed) {cover_note}", fontsize=14)

    for a in ax:
        a.axis("off")
    fig.savefig(path, dpi=60, bbox_inches="tight")
    plt.close(fig)


def save_field(outdir, filipin_path, masks, res, img=None):
    """Write one field's mask + overlay + per-cell CSV (with line/treatment/field/
    image columns) into the canonical layout. Does NOT rebuild the aggregate."""
    for sub in ("percell", "masks", "overlays"):
        os.makedirs(os.path.join(outdir, sub), exist_ok=True)
    stem = os.path.splitext(os.path.basename(filipin_path))[0]
    line, treat, field = parse_meta(stem)
    tifffile.imwrite(os.path.join(outdir, "masks", f"{stem}.masks.tif"),
                     masks.astype(np.uint16))
    df = res["df"].copy()
    df.insert(0, "image", stem)
    df.insert(0, "field", field)
    df.insert(0, "treatment", treat)
    df.insert(0, "line", line)
    df.to_csv(os.path.join(outdir, "percell", f"{stem}.csv"), index=False)
    if img is None:
        from segment import load_gray
        img = load_gray(filipin_path)
    cover = f"{100*res['df'].puncta_cover_fixed.mean():.1f}%"
    save_overlay(os.path.join(outdir, "overlays", f"{stem}.png"),
                 img, masks, res, cover)
    return stem


def rebuild_aggregate(outdir):
    """Rebuild all_cells.csv + summaries from every percell/*.csv. Returns #fields."""
    files = sorted(glob.glob(os.path.join(outdir, "percell", "*.csv")))
    if not files:
        return 0
    master = pd.concat([pd.read_csv(f) for f in files], ignore_index=True)
    master.to_csv(os.path.join(outdir, "all_cells.csv"), index=False)
    metrics = [c for c in _METRICS if c in master.columns]
    per_field = (master.groupby(["line", "treatment", "field"])
                 .agg(n_cells=("cell_id", "size"),
                      **{c: (c, "mean") for c in metrics}).reset_index())
    per_field.to_csv(os.path.join(outdir, "summary_per_field.csv"), index=False)
    per_cond = (per_field.groupby(["line", "treatment"])
                .agg(n_fields=("field", "size"), n_cells=("n_cells", "sum"),
                     **{c: (c, "mean") for c in metrics}).round(3))
    per_cond.to_csv(os.path.join(outdir, "summary_per_condition.csv"))
    return len(files)
