# Filipin Quant

Per-cell quantification of free cholesterol in **Filipin**-stained fluorescence
images (Opera Phenix, raw 16-bit TIFF). For every cell, Filipin Quant measures
how much signal sits at the **plasma membrane**, how much is in the
**cytoplasm**, and how much is packed into bright **puncta**.

It comes in two forms that share one analysis engine, so both give identical
numbers for the same settings:

- **Interactive GUI** (napari): load an image, auto-segment the cells, fix
  mistakes by hand, measure, save.
- **Batch script**: process a whole plate in one command.

---

## Contents

- [How it works](#how-it-works)
- [Installation](#installation)
- [Quick start](#quick-start)
- [Using the GUI](#using-the-gui)
- [Parameters](#parameters)
- [Output files](#output-files)
- [Input data](#input-data)
- [Troubleshooting](#troubleshooting)
- [Repository layout](#repository-layout)

---

## How it works

```
Filipin TIFF ──► Segment cells ──► (optional) hand-correct ──► Membrane / cytoplasm / puncta ──► Per-cell CSV
                 Cellpose-SAM          napari paint tools          per-cell measurements          + summaries
```

1. **Segmentation.** Cells are found on the Filipin channel itself with
   [Cellpose-SAM](https://github.com/MouseLand/cellpose) (Cellpose 4). Small
   fragments are merged into their neighbours, and isolated debris is removed.
2. **Membrane.** Each cell's outer ring is treated as membrane. A stricter
   "bright membrane" subset keeps only ring pixels clearly brighter than the
   image's typical cytoplasm level.
3. **Puncta.** After background removal (top-hat), cytoplasm pixels above a
   fixed intensity threshold count as puncta. A halo next to the membrane is
   excluded so the bright cortex is never counted as puncta.
4. **Metrics.** Mean intensities, membrane/cytoplasm ratio, % area covered by
   bright membrane and by puncta, and each compartment's share of the
   whole-cell signal.

---

## Installation

Requires [conda](https://docs.conda.io/en/latest/miniconda.html) (Miniconda
or Anaconda). An NVIDIA GPU is strongly recommended: segmentation takes about
15 s per image on GPU versus several minutes on CPU.

```bash
git clone https://github.com/Garvan-Data-Science-Platform/filipin-quant.git
cd filipin-quant
conda env create -f environment.yml
conda activate cell_quant_sam
```

> **GPU / CUDA note.** `environment.yml` pins the CUDA 12.8 (`cu128`) build of
> PyTorch. If your GPU driver is older, edit the `--extra-index-url` line to
> match your CUDA version (for example `cu121`) before creating the
> environment. Run `nvidia-smi` to see which CUDA version your driver supports.

Check the install (should print `torch cuda: True` on a GPU machine):

```bash
python -c "import napari, cellpose, torch; print('torch cuda:', torch.cuda.is_available())"
```

The first segmentation downloads the Cellpose-SAM model weights (a few hundred
MB), so the first run needs internet access and may look slow.

**Windows users:** see [`windows/README_WINDOWS.md`](windows/README_WINDOWS.md)
for a double-click desktop launcher, including how to install for a
non-admin user.

---

## Quick start

**Interactive (one image at a time):**

```bash
conda activate cell_quant_sam
python filipin_gui.py
```

**Whole plate (batch):** put your images under `data/` (one sub-folder per
plate or cell line, see [Input data](#input-data)), then:

```bash
conda activate cell_quant_sam
python batch_all.py          # writes everything to results_all/
```

---

## Using the GUI

A napari window opens with the **Filipin quant** control panel on the right.
Work down the panel from top to bottom:

1. **Load image.** Click **Browse**, pick a `...-ch1sk1fk1fl1.tif` file (the
   Filipin channel), then **Load image**.
2. **Segment.** Click **Segment**. Each detected cell gets its own colour.
   The defaults usually work (see [Parameters](#parameters)).
3. **Hand-correct.** Select the `cells` layer and use napari's paint, fill and
   erase tools (top-left toolbar) to add missed cells or delete bad masks.
   This is the only manual step, and every later number depends on it.
4. **Detect membrane.** Draws a **green** ring on each cell edge and marks the
   genuinely bright membrane in **magenta**.
5. **Detect puncta.** Marks bright cytoplasmic spots in **cyan**.
6. **Compute full metrics.** Builds the complete per-cell table.
7. **Save / update results.** Choose a results folder (for example
   `results_all`) and save. This writes the mask, per-cell CSV and QC overlay
   for this image and rebuilds the combined tables.

**Tips**

- After editing masks or changing a parameter, just click **Detect** again; it
  re-runs in a second or two.
- **Batch process folder** (bottom box) runs every Filipin TIFF in a folder
  with the current panel settings.
- **Refine after batch:** run a batch, then open any problem image, fix it and
  **Save**. Only that image's rows are replaced in the combined tables.

**What a good fit looks like:** green rings are continuous around each cell,
cyan spots sit inside cells rather than tracing the edges, and every visible
cell has a coloured mask.

---

## Parameters

The defaults were tuned on the original plates and usually need no change.

> **Keep settings fixed across an experiment.** Once you settle on the
> segmentation and puncta threshold values, use the same values for every
> image. Compare each treatment with its own vehicle control within the same
> cell line.

### Segmentation

| Parameter | Default | When to change it |
|---|---|---|
| diameter (0 = auto) | 100 | Expected cell size in pixels. Raise it if large cells are missed or split; lower it if cells merge. Auto is unreliable for these images. |
| flow_threshold | 0.90 | Higher keeps **more** cells. Raise toward 1.0 if faint cells are missed. |
| cellprob_threshold | -1.0 | Lower (toward -6) to grow masks and catch faint cells; raise (toward +6) to shrink masks and drop specks. |
| merge fragments < frac | 0.35 | Pieces smaller than this fraction of the median cell area are merged into a neighbour. |
| remove isolated tiny cells | on | Removes small lone debris. Turn off for genuinely sparse cells. |
| use GPU | on | Turn off only if the GPU is busy or unavailable. |

### Membrane

| Parameter | Default | When to change it |
|---|---|---|
| membrane ring width (px) | 6 | Thickness of the outer ring treated as membrane. |
| bright memb floor k | 3.0 | How much brighter than the image's cytoplasm level (in MADs) a ring pixel must be to count as bright membrane. Higher is stricter. |
| bright memb band (px) | 3 | How close to the cell border the bright-membrane search stays. Rarely needs changing. |

### Puncta

| Parameter | Default | When to change it |
|---|---|---|
| tophat radius (px) | 15 | Roughly the size of the largest puncta. Increase for larger accumulations. |
| opening radius (px) | 2 | Raise slightly if thin membrane signal leaks into the puncta mask. |
| membrane halo excl. (px) | 4 | Buffer next to the membrane excluded from puncta. Raise if cyan hugs the cell edge. |
| fixed threshold | 150 | **The key puncta setting.** Background-subtracted intensity above which a pixel counts as puncta. Lower gives more coverage, higher gives less. |

---

## Output files

Everything is written to one results folder (`results_all/` for the batch
script):

| File | Contents |
|---|---|
| `percell/<image>.csv` | Per-cell metrics for one image. **This is the source of truth.** |
| `all_cells.csv` | Every cell from every image, combined. |
| `summary_per_field.csv` | One row per image (averages). |
| `summary_per_condition.csv` | One row per cell line x treatment. |
| `masks/<image>.masks.tif` | Label mask (one integer ID per cell). |
| `overlays/<image>.png` | QC overlay image. |
| `params.json` | Exact parameter values used (GUI saves). |

The combined files are always rebuilt from `percell/`, so re-saving one
corrected image updates the totals without touching anything else.

**Key per-cell columns**

| Column | Meaning |
|---|---|
| `membrane_mean` | Mean Filipin intensity in the membrane ring. |
| `cytoplasm_mean` | Mean intensity in the cytoplasm (puncta excluded). |
| `membrane_to_cyto_ratio` | `membrane_mean / cytoplasm_mean`. |
| `bright_membrane_cover` | Fraction of the ring that is bright membrane. |
| `puncta_cover_fixed` | Fraction of cytoplasm covered by puncta (main puncta readout). |
| `membrane_frac`, `cytoplasm_frac`, `puncta_frac` | Each compartment's share of the whole-cell signal. |

---

## Input data

- **Format:** raw **16-bit single-channel TIFFs** exported from Harmony (Opera
  Phenix). Do not use 8-bit or contrast-stretched exports; their intensities
  are not comparable between images.
- **Channels:** `ch1` = Filipin (quantified), `ch2` = brightfield (not used).
  Only files ending in `-ch1sk1fk1fl1.tif` (or `.tiff`) are processed.
- **Acquisition:** keep exposure, gain and imaging order identical across all
  wells. Filipin photobleaches quickly.
- **Folder layout for `batch_all.py`:** `data/<any sub-folder>/<image files>`.
- **Condition labels:** the cell line, treatment and field are read from the
  file name by `parse_meta()` in `results_io.py`. If your naming differs,
  edit that function so the summary tables group correctly.

Image data and results are not part of this repository.

---

## Troubleshooting

| Symptom | Fix |
|---|---|
| Segmentation "Killed" or very slow | GPU is busy or not detected. Check `torch.cuda.is_available()`, or untick **use GPU**. |
| Large cells missed or split | Raise **diameter**. |
| Faint cells missed | Raise **flow_threshold** (toward 1.0) or lower **cellprob_threshold**. |
| Cyan puncta along the membrane | Raise **membrane halo excl.** or **opening radius**. |
| Too much or too little puncta | Adjust **fixed threshold** (and keep it the same for the whole experiment). |
| First run seems to hang | It is downloading the Cellpose-SAM weights; wait for it to finish. |

---

## Repository layout

| File | Purpose |
|---|---|
| `filipin_gui.py` | napari GUI (single image + batch folder). |
| `batch_all.py` | Command-line batch over everything in `data/`. |
| `segment.py` | Cellpose-SAM segmentation, fragment merging, debris removal. |
| `filipin_quant.py` | Membrane / cytoplasm / puncta measurements (`analyze()`). |
| `results_io.py` | Writes per-image results and rebuilds combined tables. |
| `render_overlays.py` | Stand-alone QC overlay rendering from a saved mask. |
| `environment.yml` | Pinned conda environment (`cell_quant_sam`). |
| `windows/` | Windows desktop launcher and shortcut scripts. |
