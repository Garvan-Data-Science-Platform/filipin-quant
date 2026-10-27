#!/usr/bin/env python3
"""
Render QC overlays + save membrane/cytoplasm/puncta masks for a segmented field.

Reuses the exact measurement logic from filipin_quant.py so what you SEE is
what was measured. Reads the *.masks.tif written by --save-masks.

Outputs (prefix = --out stem):
  <p>_overlay.png        4-panel QC figure
  <p>_membrane.tif       binary shell-membrane mask (uint8 0/255)
  <p>_cytoplasm.tif      binary diffuse-cytoplasm mask, puncta excluded (uint8 0/255)
  <p>_puncta.tif         binary puncta mask, fixed-threshold coverage (uint8 0/255)

Usage:
  python render_overlays.py --filipin field-ch1.tif --masks field.masks.tif --out field
"""
import argparse

import numpy as np
import tifffile
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from skimage.segmentation import find_boundaries
from skimage.morphology import binary_dilation, disk

from filipin_quant import load_gray, analyze


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--filipin", required=True)
    ap.add_argument("--masks", required=True, help="*.masks.tif from filipin_quant.py")
    ap.add_argument("--out", required=True, help="output prefix")
    ap.add_argument("--ring", type=int, default=6)
    ap.add_argument("--tophat", type=int, default=15)
    ap.add_argument("--punc-open", type=int, default=2)
    ap.add_argument("--memexclude", type=int, default=4)
    ap.add_argument("--punc-thresh", type=float, default=150.0)
    ap.add_argument("--punc-mode", choices=["adaptive", "fixed"], default="fixed",
                    help="Which puncta mask to display in the overlay")
    args = ap.parse_args()

    filipin = load_gray(args.filipin)
    masks = tifffile.imread(args.masks).astype(np.int32)

    res = analyze(filipin, masks, ring_px=args.ring,
                  tophat_radius=args.tophat, membrane_exclude_px=args.memexclude,
                  punc_open=args.punc_open, punc_fixed_thresh=args.punc_thresh)
    membrane, cyto = res["membrane"], res["cytoplasm"]
    pmask = res["pmask_adaptive"] if args.punc_mode == "adaptive" else res["pmask"]
    df = res["df"]
    n_cells = int(masks.max())
    cover = 100 * (df.puncta_cover_adaptive.mean() if args.punc_mode == "adaptive"
                   else df.puncta_cover_fixed.mean())
    print(f"{args.out}: {n_cells} cells, mean puncta cover ({args.punc_mode}) {cover:.1f}%")

    # --- save masks ---
    tifffile.imwrite(f"{args.out}_membrane.tif", (membrane * 255).astype(np.uint8))
    tifffile.imwrite(f"{args.out}_cytoplasm.tif", (cyto * 255).astype(np.uint8))
    tifffile.imwrite(f"{args.out}_puncta.tif", (pmask * 255).astype(np.uint8))

    # --- 4-panel QC figure ---
    vmax = np.percentile(filipin, 99.5)
    base = np.clip(filipin / vmax, 0, 1)
    fig, ax = plt.subplots(2, 2, figsize=(16, 13))
    for a in ax.ravel():
        a.axis("off")

    ax[0, 0].imshow(filipin, cmap="gray", vmax=vmax)
    ax[0, 0].set_title("Filipin (raw display)")

    ov = np.dstack([base] * 3)
    ov[find_boundaries(masks, mode="outer")] = [1, 1, 0]  # yellow cell borders
    ax[0, 1].imshow(ov)
    ax[0, 1].set_title(f"Cell segmentation ({n_cells} cells)")

    seg = np.zeros((*masks.shape, 3))
    seg[..., 0][cyto] = 0.35        # cytoplasm -> dim red
    seg[membrane] = [0, 0.9, 0]     # shell membrane -> green
    ax[1, 0].imshow(base, cmap="gray")
    ax[1, 0].imshow(seg, alpha=0.55)
    ax[1, 0].set_title("Membrane (green) / cytoplasm (red)")

    pov = np.dstack([base] * 3)
    pov[binary_dilation(pmask, disk(1))] = [0, 1, 1]  # puncta -> cyan
    ax[1, 1].imshow(pov)
    ax[1, 1].set_title(f"Puncta coverage ({args.punc_mode}): "
                       f"mean {cover:.1f}% of cytoplasm")

    fig.tight_layout()
    fig.savefig(f"{args.out}_overlay.png", dpi=120, bbox_inches="tight")
    print(f"  wrote {args.out}_overlay.png + membrane/cytoplasm/puncta .tif")


if __name__ == "__main__":
    main()
