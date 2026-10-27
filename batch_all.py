#!/usr/bin/env python3
"""Quantify EVERY Filipin field in data/ with the finalised pipeline:
   Cellpose-SAM (diam 100, flow 0.9, cellprob -1.0) + fragment-merge +
   largest-component + isolated-cell removal  ->  analyze() metrics.

Run in the cell_quant_sam (cellpose>=4) env.

Writes the canonical result set (see results_io) under results_all/:
  percell/<stem>.csv (source of truth), masks/, overlays/,
  all_cells.csv, summary_per_field.csv, summary_per_condition.csv.

Because per-field CSVs are the source of truth, refining a single field later
in the GUI (Save) overwrites that field and rebuilds the aggregate — no need to
re-run the whole plate.
"""
import glob
import os

from segment import load_gray, segment_sam, merge_fragments, remove_isolated
from filipin_quant import analyze
from results_io import save_field, rebuild_aggregate

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "data")
OUT = os.path.join(HERE, "results_all")


def main():
    files = sorted(set(
        glob.glob(f"{DATA}/*/*-ch1sk1fk1fl1.tif")
        + glob.glob(f"{DATA}/*/*-ch1sk1fk1fl1.tiff")
        + glob.glob(f"{DATA}/R172H*/*-ch1sk1fk1fl1.tif")
        + glob.glob(f"{DATA}/R172H*/*-ch1sk1fk1fl1.tiff")))
    print(f"{len(files)} Filipin fields", flush=True)
    for i, f in enumerate(files):
        img = load_gray(f)
        masks = segment_sam(img, diameter=100)         # flow 0.9 / cellprob -1.0 defaults
        masks, _ = merge_fragments(masks)              # + keep_largest_component inside
        masks, _ = remove_isolated(masks)
        res = analyze(img, masks)
        stem = save_field(OUT, f, masks, res, img=img)
        df = res["df"]
        print(f"[{i+1}/{len(files)}] {stem}: {len(df)} cells, "
              f"membrane {df.membrane_mean.mean():.0f}, "
              f"bright-memb {100*df.bright_membrane_cover.mean():.0f}%, "
              f"puncta {100*df.puncta_cover_fixed.mean():.1f}%", flush=True)

    n = rebuild_aggregate(OUT)
    import pandas as pd
    cond = pd.read_csv(f"{OUT}/summary_per_condition.csv")
    print(f"\n===== SUMMARY (within-line), {n} fields =====")
    print(cond.to_string(index=False))
    print("DONE")


if __name__ == "__main__":
    main()
