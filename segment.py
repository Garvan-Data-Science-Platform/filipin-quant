#!/usr/bin/env python3
"""Segment the Filipin channel with Cellpose-SAM (v4) and merge small fragments.

Chosen backend after benchmarking: SAM with an explicit diameter (default 70)
gives the best masks across all three cell lines at ONE setting (no per-line
tuning). SAM occasionally leaves thin fragments between cells; `merge_fragments`
absorbs each small piece into the neighbour it shares the most border with.

Run in the cell_quant_sam (cellpose>=4) env. Writes <out>.masks.tif (+ overlay).

Usage:
  conda run -n cell_quant_sam python segment.py --filipin field-ch1.tif --out field
"""
import argparse

import numpy as np
import tifffile
from scipy import ndimage as ndi
from skimage.segmentation import relabel_sequential


def load_gray(path):
    img = tifffile.imread(path)
    if img.ndim == 3 and img.shape[-1] in (3, 4):
        img = img[..., :3].mean(axis=-1)
    return img.astype(np.float32)


def segment_sam(img, diameter=100, flow_threshold=0.9, cellprob_threshold=-1.0,
                gpu=True):
    from cellpose import models
    model = models.CellposeModel(gpu=gpu)          # v4 single generalist model
    # NB: eval normalizes the input array IN PLACE — pass a copy so the caller's
    # raw image is preserved for downstream intensity quantification.
    masks = model.eval(img.copy(), diameter=diameter, flow_threshold=flow_threshold,
                       cellprob_threshold=cellprob_threshold)[0]
    return masks.astype(np.int32)


def merge_fragments(masks, min_area_frac=0.35, drop_abs_area=150):
    """Merge small fragments into the adjacent cell with the longest shared border.

    A fragment is a label whose area < `min_area_frac` * median cell area.
    If a fragment has no cell neighbour and is tiny (< `drop_abs_area`), it is
    treated as debris and removed; otherwise the smallest neighbours are kept.
    """
    masks = masks.copy()
    labels, counts = np.unique(masks[masks > 0], return_counts=True)
    if labels.size == 0:
        return masks, 0
    thresh = min_area_frac * np.median(counts)
    small = [int(l) for l, c in sorted(zip(labels, counts), key=lambda x: x[1])
             if c < thresh]

    slices = ndi.find_objects(masks)               # index by (label - 1)
    merged = 0
    for lab in small:
        sl = slices[lab - 1] if lab - 1 < len(slices) else None
        if sl is None:
            continue
        pad = tuple(slice(max(0, s.start - 3), min(d, s.stop + 3))
                    for s, d in zip(sl, masks.shape))
        sub = masks[pad]                            # view into masks
        frag = sub == lab
        if not frag.any():
            continue
        ring = ndi.binary_dilation(frag, iterations=2) & ~frag
        neigh = sub[ring]
        neigh = neigh[(neigh > 0) & (neigh != lab)]
        if neigh.size == 0:
            if frag.sum() < drop_abs_area:          # isolated debris -> drop
                sub[frag] = 0
                merged += 1
            continue
        vals, cnts = np.unique(neigh, return_counts=True)
        sub[frag] = int(vals[np.argmax(cnts)])      # longest shared border wins
        merged += 1

    masks = keep_largest_component(masks)           # one connected blob per label
    masks, _, _ = relabel_sequential(masks)
    return masks, merged


def keep_largest_component(masks):
    """Ensure every label is a single connected blob: if a label's mask is split
    into several pieces (a segmentation/merge artefact that also breaks the
    membrane ring), keep only its largest component and drop the rest.
    """
    out = masks.copy()
    full = np.ones((3, 3), int)                     # 8-connectivity
    for idx, sl in enumerate(ndi.find_objects(masks)):
        lab = idx + 1
        if sl is None:
            continue
        sub = out[sl]
        cell = sub == lab
        comp, n = ndi.label(cell, structure=full)
        if n <= 1:
            continue
        sizes = np.bincount(comp.ravel())
        sizes[0] = 0
        largest = int(sizes.argmax())
        sub[cell & (comp != largest)] = 0           # drop the smaller pieces
    return out


def remove_isolated(masks, max_area_frac=0.4, gap_px=6):
    """Remove small cells that have NO neighbouring cell within `gap_px`.

    Targets isolated debris/artefacts (a tiny object sitting alone on background)
    while keeping genuine cells: a label is removed only if BOTH its area is
    < `max_area_frac` * median cell area AND nothing else lies within `gap_px`.
    Large isolated cells (e.g. sparse fields) are kept.
    """
    labels, counts = np.unique(masks[masks > 0], return_counts=True)
    if labels.size == 0:
        return masks, 0
    thresh = max_area_frac * np.median(counts)
    slices = ndi.find_objects(masks)
    removed = 0
    for lab, cnt in zip(labels, counts):
        if cnt >= thresh:
            continue
        sl = slices[lab - 1]
        if sl is None:
            continue
        pad = tuple(slice(max(0, s.start - gap_px - 1), min(d, s.stop + gap_px + 1))
                    for s, d in zip(sl, masks.shape))
        sub = masks[pad]
        frag = sub == lab
        ring = ndi.binary_dilation(frag, iterations=gap_px) & ~frag
        neigh = sub[ring]
        if not np.any((neigh > 0) & (neigh != lab)):     # no neighbour -> isolated
            sub[frag] = 0
            removed += 1
    if removed:
        masks, _, _ = relabel_sequential(masks)
    return masks, removed


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--filipin", required=True)
    ap.add_argument("--out", required=True, help="output prefix")
    ap.add_argument("--diam", type=float, default=100)
    ap.add_argument("--flow", type=float, default=0.9)
    ap.add_argument("--cellprob", type=float, default=-1.0)
    ap.add_argument("--min-area-frac", type=float, default=0.35)
    ap.add_argument("--no-merge", action="store_true")
    args = ap.parse_args()

    img = load_gray(args.filipin)
    masks = segment_sam(img, diameter=args.diam, flow_threshold=args.flow,
                        cellprob_threshold=args.cellprob)
    n0 = int(masks.max())
    if not args.no_merge:
        masks, merged = merge_fragments(masks, min_area_frac=args.min_area_frac)
        masks, iso = remove_isolated(masks)
        print(f"{args.out}: {n0} raw -> {int(masks.max())} cells "
              f"({merged} fragments merged, {iso} isolated removed)")
    else:
        print(f"{args.out}: {n0} cells (no merge)")

    tifffile.imwrite(f"{args.out}.masks.tif", masks.astype(np.uint16))


if __name__ == "__main__":
    main()
