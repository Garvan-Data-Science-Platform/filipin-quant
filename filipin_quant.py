#!/usr/bin/env python3
"""
Filipin (free-cholesterol) per-cell quantification for Opera Phenix data.

Measures, per segmented cell:
  - membrane intensity   (bright cortical ring)
  - cytoplasm intensity  (interior, with puncta excluded)
  - puncta count         (LoG blob detection)

IMPORTANT
---------
Run this on the RAW 16-bit single-channel TIFFs exported from Harmony,
NOT on 8-bit RGB PNG display exports. The PNGs are contrast-stretched and
are not quantitatively meaningful.

Segmentation is done on the FILIPIN channel (membrane signal gives clean
cell borders), which works better than brightfield for confluent monolayers.
Pass a brightfield image too if you want to try it instead.

Deps:
    pip install cellpose scikit-image tifffile numpy pandas scipy

Usage:
    python filipin_quant.py --filipin path/to/filipin.tif --out results.csv
    python filipin_quant.py --filipin f.tif --diam 60 --out results.csv
"""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import tifffile
from scipy import ndimage as ndi
from skimage.feature import blob_log
from skimage.morphology import (
    binary_erosion, disk, opening, remove_small_objects, white_tophat,
)
from skimage.segmentation import find_boundaries
from skimage.measure import regionprops


def load_gray(path):
    """Load a single-channel image; if RGB sneaks in, collapse to luminance."""
    img = tifffile.imread(path) if str(path).lower().endswith((".tif", ".tiff")) \
        else __import__("skimage.io", fromlist=["imread"]).imread(path)
    if img.ndim == 3 and img.shape[-1] in (3, 4):
        img = img[..., :3].mean(axis=-1)
    return img.astype(np.float32)


def subtract_background(img, radius=50, downsample=4):
    """Rolling-ball-like background removal via a large morphological opening.

    white_tophat(img, disk(radius)) = img - opening(img, disk(radius)). A full-res
    disk(50) opening on a 2160^2 image costs ~60 s. The background (the opening) is
    smooth, so we estimate it on a `downsample`x-decimated image with a proportionally
    smaller disk and upsample — ~downsample^2 faster (~4 s) with negligible change to
    the subtracted result. Set downsample=1 for the exact full-res tophat.
    """
    img = img.astype(np.float32)
    if downsample and downsample > 1:
        from skimage.transform import resize
        small = img[::downsample, ::downsample]
        r = max(1, int(round(radius / downsample)))
        bg = resize(opening(small, disk(r)), img.shape, order=1,
                    preserve_range=True, anti_aliasing=False)
        return np.clip(img - bg, 0, None).astype(np.float32)
    return white_tophat(img, footprint=disk(radius)).astype(np.float32)


def segment_cells(img, diameter, use_gpu=False, flow_threshold=0.4,
                  cellprob_threshold=0.0):
    """Whole-cell segmentation on the Filipin (membrane) channel with Cellpose 3.

    Uses CellposeDenoiseModel (cyto3 + integrated denoising), which suits the
    noisy Filipin signal. `diameter` controls how confluent cells get split;
    pass an explicit value if auto-estimation over-merges the monolayer.

    For the dense positive-control monolayer, diam=50 flow_threshold=0.9 gives
    full coverage without over-merging. The cleaner syngeneic-43 ring works well
    at the diam=60 / flow=0.4 defaults.
    """
    from cellpose import denoise
    model = denoise.CellposeDenoiseModel(
        gpu=use_gpu, model_type="cyto3", restore_type="denoise_cyto3",
    )
    # Single grayscale channel -> channels=[0, 0]. eval returns
    # (masks, flows, styles, imgs_denoised) in the denoise model.
    masks, _flows, _styles, _denoised = model.eval(
        img, diameter=diameter, channels=[0, 0],
        flow_threshold=flow_threshold, cellprob_threshold=cellprob_threshold,
    )
    return masks.astype(np.int32)


def detect_puncta(img, mask=None, min_sigma=2, max_sigma=6, num_sigma=5,
                  threshold_rel=0.1):
    """LoG blob detection on a tophat-enhanced image, restricted to `mask` if given.

    `mask` should be the cytoplasm interior (membrane excluded) so the detector
    does not fire on the bright plasma-membrane ridge. Sigmas are in pixels and
    should scale with image resolution (defaults suit full-res 2160px fields).
    """
    enhanced = white_tophat(img, footprint=disk(int(max_sigma * 2)))
    if mask is not None:
        enhanced = enhanced * mask  # suppress everything outside the search region
    if enhanced.max() > 0:
        enhanced = enhanced / enhanced.max()
    blobs = blob_log(enhanced, min_sigma=min_sigma, max_sigma=max_sigma,
                     num_sigma=num_sigma, threshold=None, threshold_rel=threshold_rel)
    return blobs[:, :2].astype(int) if len(blobs) else np.empty((0, 2), int)


def puncta_mask_from_coords(coords, shape, radius=2):
    m = np.zeros(shape, bool)
    for r, c in coords:
        if 0 <= r < shape[0] and 0 <= c < shape[1]:
            m[r, c] = True
    return ndi.binary_dilation(m, structure=disk(radius))


def cell_regions(cell, ring_px):
    """Return (membrane, interior) for ONE cell.

    membrane : outer `ring_px` morphological ring of the mask — a complete, closed
               annulus tracing the whole cell boundary. Because segmentation is done
               on the Filipin signal, the mask border sits on the membrane.
    interior : the eroded interior (outer ring removed) = cytoplasm territory.

    The erosion is done on a crop PADDED with background, otherwise binary_erosion
    treats the bounding-box edge as interior and the ring breaks into arcs wherever
    the cell touches its bbox (visible as gaps on colony edges / isolated cells).
    """
    pad = ring_px + 1
    cp = np.pad(cell, pad, mode="constant", constant_values=False)
    interior = binary_erosion(cp, disk(ring_px))[pad:-pad, pad:-pad]
    membrane = cell & ~interior
    return membrane, interior


def compute_bg(filipin_raw):
    """Precompute the expensive, parameter-independent rolling-ball background so
    an interactive GUI can cache it and re-run detection quickly on threshold changes."""
    return subtract_background(filipin_raw)


def analyze(filipin_raw, masks, ring_px=6, tophat_radius=15,
            membrane_exclude_px=4, punc_open=2, punc_fixed_thresh=150.0,
            punc_adapt_k=2.0, bright_floor_k=3.0, bright_band_px=3,
            bright_min_arc=15, bg=None):
    """Full per-cell analysis on the Filipin channel.

    Membrane = the outer `ring_px` MORPHOLOGICAL RING of each cell mask. Because
    segmentation is done on the membranous Filipin signal, the mask border sits on
    the membrane, so the ring gives complete, gap-free coverage.

    Puncta (intracellular free-cholesterol accumulations) are summarised as
    %-AREA COVERAGE of the cytoplasm, not blob counts, because their size and
    intensity vary widely. A FIXED absolute threshold on the enhanced signal is the
    primary metric (comparable across images given identical acquisition); a per-cell
    ADAPTIVE threshold is also reported for reference. Puncta and membrane are
    separated by SHAPE: a morphological OPENING (`punc_open`) removes thin membrane
    residue while keeping compact blobs, and a thin `membrane_exclude_px` halo around
    the ring is excluded from the cytoplasm.

    Returns a dict: per-cell DataFrame + global masks for QC overlays.
    Per-cell work uses bounding-box crops (ndi.find_objects) for speed.
    """
    if bg is None:
        bg = subtract_background(filipin_raw)            # rolling-ball background
    enhanced = white_tophat(bg, disk(tophat_radius))     # local bright features
    if punc_open > 0:
        enhanced = opening(enhanced, disk(punc_open))    # drop thin membrane, keep blobs

    slices = ndi.find_objects(masks)
    memb_full = np.zeros(masks.shape, bool)
    interior_full = np.zeros(masks.shape, bool)
    cache = {}
    for idx, sl in enumerate(slices):
        lab = idx + 1
        if sl is None:
            continue
        cellsub = masks[sl] == lab
        mm, interior = cell_regions(cellsub, ring_px)
        memb_full[sl] |= mm
        interior_full[sl] |= interior
        cache[lab] = (sl, mm, cellsub)

    # BRIGHT ('real') membrane: pixels in the ring that are genuinely bright
    # (>= a global cytoplasm-derived floor), restricted to a band around the
    # cell-cell / cell-background interface so bright juxta-membrane puncta
    # sitting inside the ring are dropped; then connected-arc filtered (arcs
    # allowed — the membrane need not be a closed loop). Complementary to the
    # geometry-ring intensity metrics, which are left untouched.
    cyto_vals = bg[interior_full]
    if cyto_vals.size:
        med = float(np.median(cyto_vals))
        mad = float(np.median(np.abs(cyto_vals - med)))
        bright_floor = med + bright_floor_k * 1.4826 * (mad if mad > 0 else cyto_vals.std())
    else:
        bright_floor = np.inf
    boundary_band = ndi.binary_dilation(find_boundaries(masks, mode="outer"),
                                        structure=disk(bright_band_px))
    bright_full = memb_full & (bg >= bright_floor) & boundary_band
    bright_full = remove_small_objects(bright_full, bright_min_arc, connectivity=2)

    # cytoplasm territory = interior minus a halo around the membrane ring
    cyto_terr_full = interior_full & ~ndi.binary_dilation(
        memb_full, structure=disk(membrane_exclude_px))
    punc_fixed_full = cyto_terr_full & (enhanced >= punc_fixed_thresh)

    cyto_diffuse_full = np.zeros(masks.shape, bool)
    punc_adapt_full = np.zeros(masks.shape, bool)
    rows = []
    for lab, (sl, mm, cellsub) in cache.items():
        bgsub = bg[sl]
        ensub = enhanced[sl]
        terr = cyto_terr_full[sl] & cellsub               # this cell's cytoplasm territory
        p_fixed = terr & (ensub >= punc_fixed_thresh)
        vals = ensub[terr]
        if vals.size:
            p_adapt = terr & (ensub >= vals.mean() + punc_adapt_k * vals.std())
        else:
            p_adapt = np.zeros_like(terr)
        punc_adapt_full[sl] |= p_adapt
        cyto_diffuse = terr & ~p_fixed                    # diffuse cytoplasm (puncta removed)
        cyto_diffuse_full[sl] |= cyto_diffuse

        terr_area = int(terr.sum())
        ring_area = int(mm.sum())
        bright_cell = bright_full[sl] & mm                # bright 'real' membrane in this ring
        mm_mean = float(bgsub[mm].mean()) if mm.any() else np.nan
        cyto_mean = float(bgsub[cyto_diffuse].mean()) if cyto_diffuse.any() else np.nan
        # integrated signal per compartment -> fractions of whole-cell Filipin.
        # These are dimensionless (0-1), so they are directly comparable across
        # images, cell lines and treatments regardless of absolute staining level.
        whole_int = float(bgsub[cellsub].sum())
        memb_int = float(bgsub[mm].sum())
        cyto_int = float(bgsub[cyto_diffuse].sum())
        punc_int = float(bgsub[p_fixed].sum())
        rows.append({
            "cell_id": lab,
            "area_px": int(cellsub.sum()),
            "membrane_area_px": int(mm.sum()),
            "cyto_area_px": terr_area,
            # --- absolute intensities (comparable given identical acquisition) ---
            "membrane_mean": mm_mean,             # morphological ring (all ring px)
            "bright_membrane_mean": (float(bgsub[bright_cell].mean())
                                     if bright_cell.any() else np.nan),
            "cytoplasm_mean": cyto_mean,
            "whole_cell_mean": float(bgsub[cellsub].mean()),
            # --- normalized / dimensionless (comparable across everything) ---
            "membrane_to_cyto_ratio": (mm_mean / cyto_mean
                                       if cyto_mean and cyto_mean > 0 else np.nan),
            "bright_membrane_cover": (bright_cell.sum() / ring_area
                                      if ring_area else np.nan),
            "membrane_frac": (memb_int / whole_int if whole_int > 0 else np.nan),
            "cytoplasm_frac": (cyto_int / whole_int if whole_int > 0 else np.nan),
            "puncta_frac": (punc_int / whole_int if whole_int > 0 else np.nan),
            "puncta_cover_fixed": (p_fixed.sum() / terr_area if terr_area else np.nan),
            "puncta_cover_adaptive": (p_adapt.sum() / terr_area if terr_area else np.nan),
            "puncta_intensity_sum": punc_int,
        })
    return {
        "df": pd.DataFrame(rows),
        "membrane": memb_full,                # morphological ring membrane for overlay
        "bright_membrane": bright_full,       # bright 'real' membrane (boundary-band)
        "cytoplasm": cyto_diffuse_full,
        "pmask": punc_fixed_full,             # puncta overlay = fixed-threshold coverage
        "pmask_adaptive": punc_adapt_full,
    }


def quantify(filipin_raw, masks, **kw):
    """Thin wrapper: returns the per-cell DataFrame."""
    return analyze(filipin_raw, masks, **kw)["df"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--filipin", required=True, help="Raw 16-bit Filipin TIFF")
    ap.add_argument("--out", default="filipin_results.csv")
    ap.add_argument("--masks", default=None,
                    help="Precomputed label-mask TIFF (skip segmentation; e.g. from segment.py)")
    ap.add_argument("--diam", type=float, default=None,
                    help="Approx cell diameter in px (None = auto-estimate)")
    ap.add_argument("--ring", type=int, default=6, help="Membrane ring width (px)")
    ap.add_argument("--tophat", type=int, default=15,
                    help="White-tophat radius (px) for puncta enhancement")
    ap.add_argument("--punc-open", type=int, default=2,
                    help="Opening radius (px) to drop thin membrane, keep blob puncta")
    ap.add_argument("--memexclude", type=int, default=4,
                    help="Thin halo (px) around membrane excluded from cytoplasm/puncta")
    ap.add_argument("--punc-thresh", type=float, default=150.0,
                    help="Fixed absolute puncta threshold (tophat-enhanced units)")
    ap.add_argument("--punc-k", type=float, default=2.0,
                    help="Adaptive puncta threshold = cyto mean + k*SD")
    ap.add_argument("--flow", type=float, default=0.4,
                    help="Cellpose flow_threshold (raise to recover hard cells)")
    ap.add_argument("--cellprob", type=float, default=0.0,
                    help="Cellpose cellprob_threshold")
    ap.add_argument("--gpu", action="store_true")
    ap.add_argument("--save-masks", action="store_true")
    args = ap.parse_args()

    filipin = load_gray(args.filipin)
    print(f"Loaded {args.filipin}  shape={filipin.shape}  dtype range="
          f"[{filipin.min():.0f}, {filipin.max():.0f}]")
    if filipin.max() <= 255:
        print("WARNING: max <= 255 — this looks 8-bit. Use the raw 16-bit TIFF "
              "for valid intensity measurements.")

    if args.masks:
        masks = tifffile.imread(args.masks).astype(np.int32)
        print(f"Loaded masks {args.masks}: {masks.max()} cells")
    else:
        print("Segmenting on Filipin channel...")
        masks = segment_cells(filipin, args.diam, use_gpu=args.gpu,
                              flow_threshold=args.flow, cellprob_threshold=args.cellprob)
        print(f"  {masks.max()} cells found")

    df = quantify(filipin, masks, ring_px=args.ring,
                  tophat_radius=args.tophat, membrane_exclude_px=args.memexclude,
                  punc_open=args.punc_open,
                  punc_fixed_thresh=args.punc_thresh, punc_adapt_k=args.punc_k)
    df.to_csv(args.out, index=False)
    print(f"Wrote {len(df)} cells -> {args.out}")
    print(f"  mean membrane={df.membrane_mean.mean():.0f}  "
          f"cytoplasm={df.cytoplasm_mean.mean():.0f}  "
          f"membrane/cyto={df.membrane_to_cyto_ratio.mean():.2f}  "
          f"puncta cover: fixed={100*df.puncta_cover_fixed.mean():.1f}%  "
          f"adaptive={100*df.puncta_cover_adaptive.mean():.1f}%")

    if args.save_masks:
        tifffile.imwrite(Path(args.out).with_suffix(".masks.tif"),
                         masks.astype(np.uint16))


if __name__ == "__main__":
    main()
