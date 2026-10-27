#!/usr/bin/env python3
"""
Filipin quantification GUI (napari + Qt).

Single-image workflow:
  1. Load a Filipin ch1 TIFF.
  2. Segment (Cellpose-SAM + fragment-merge + isolated-cell removal).
  3. MANUALLY CORRECT the 'cells' labels layer with napari's paint/fill/erase
     tools (fix missed colonies, delete bad masks).
  4. Detect membrane / Detect puncta (each with its own overlay + readout).
  5. Compute full metrics table -> Save (mask, per-cell CSV, params).

Batch: process every ch1 TIFF in a folder with the current parameters.

Parameters for SEGMENTATION, MEMBRANE and PUNCTA are grouped in framed boxes and
resize with the dock. The engine is the shared segment.py + analyze() code, so
GUI and batch CLI give identical numbers for the same parameters.

Run (needs a display + the cell_quant_sam env):
    conda run -n cell_quant_sam python filipin_gui.py
"""
import glob
import json
import os
import warnings

warnings.filterwarnings("ignore", category=FutureWarning)   # skimage binary_erosion
warnings.filterwarnings("ignore", category=UserWarning)     # torch sparse notice

import numpy as np
import pandas as pd
import tifffile
import napari
from qtpy.QtWidgets import (
    QWidget, QVBoxLayout, QGroupBox, QFormLayout, QHBoxLayout, QPushButton,
    QDoubleSpinBox, QSpinBox, QCheckBox, QComboBox, QLineEdit, QLabel,
    QFileDialog, QScrollArea, QSizePolicy, QProgressBar, QApplication,
)
from qtpy.QtCore import Qt

from segment import segment_sam, merge_fragments, remove_isolated, load_gray
from filipin_quant import analyze, compute_bg
from results_io import save_field, rebuild_aggregate

CH1_GLOBS = ["*-ch1sk1fk1fl1.tif", "*-ch1sk1fk1fl1.tiff"]


def _dspin(val, lo, hi, step, dec=2):
    w = QDoubleSpinBox(); w.setRange(lo, hi); w.setSingleStep(step)
    w.setDecimals(dec); w.setValue(val)
    w.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)
    return w


def _ispin(val, lo, hi):
    w = QSpinBox(); w.setRange(lo, hi); w.setValue(val)
    w.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)
    return w


def _group(title):
    box = QGroupBox(title)
    box.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)
    form = QFormLayout(box)
    form.setFieldGrowthPolicy(QFormLayout.AllNonFixedFieldsGrow)  # fields fill width
    form.setLabelAlignment(Qt.AlignLeft)
    return box, form


class FilipinWidget(QWidget):
    def __init__(self, viewer):
        super().__init__()
        self.viewer = viewer
        self.image = None
        self.path = None
        self.res = None
        self._bg = None          # cached rolling-ball background (per image)

        self.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
        outer = QVBoxLayout(self)
        outer.setContentsMargins(4, 4, 4, 4)

        # ---------- Input ----------
        gb, form = _group("Input")
        self.path_edit = QLineEdit(); self.path_edit.setPlaceholderText("Filipin ch1 TIFF…")
        browse = QPushButton("Browse"); browse.clicked.connect(self._browse)
        row = QHBoxLayout(); row.addWidget(self.path_edit); row.addWidget(browse)
        form.addRow("File", self._wrap(row))
        self.btn_load = QPushButton("Load image"); self.btn_load.clicked.connect(self._load)
        form.addRow(self.btn_load)
        outer.addWidget(gb)

        # ---------- Segmentation ----------
        gb, form = _group("Segmentation")
        self.diameter = _dspin(100, 0, 500, 5, 0)
        self.flow = _dspin(0.9, 0, 3, 0.1)
        self.cellprob = _dspin(-1.0, -6, 6, 0.5)
        self.merge_frac = _dspin(0.35, 0, 1, 0.05)
        self.rm_iso = QCheckBox("remove isolated tiny cells"); self.rm_iso.setChecked(True)
        self.use_gpu = QCheckBox("use GPU"); self.use_gpu.setChecked(True)
        form.addRow("diameter (0=auto)", self.diameter)
        form.addRow("flow_threshold", self.flow)
        form.addRow("cellprob_threshold", self.cellprob)
        form.addRow("merge fragments < frac", self.merge_frac)
        form.addRow(self.rm_iso)
        form.addRow(self.use_gpu)
        self.btn_seg = QPushButton("Segment"); self.btn_seg.clicked.connect(self._segment)
        form.addRow(self.btn_seg)
        self.seg_out = self._status(); form.addRow(self.seg_out)
        outer.addWidget(gb)

        # ---------- Membrane ----------
        gb, form = _group("Membrane detection")
        self.ring = _ispin(6, 1, 50)
        self.bfloor = _dspin(3.0, 0, 10, 0.5)
        self.bband = _ispin(3, 1, 20)
        form.addRow("membrane ring width (px)", self.ring)
        form.addRow("bright memb floor k (cyto+k·MAD)", self.bfloor)
        form.addRow("bright memb band (px)", self.bband)
        self.btn_memb = QPushButton("Detect membrane"); self.btn_memb.clicked.connect(self._detect_membrane)
        form.addRow(self.btn_memb)
        self.memb_out = self._status(); form.addRow(self.memb_out)
        outer.addWidget(gb)

        # ---------- Puncta ----------
        gb, form = _group("Puncta detection")
        self.tophat = _ispin(15, 3, 60)
        self.popen = _ispin(2, 0, 10)
        self.memex = _ispin(4, 0, 40)
        self.pthresh = _dspin(150, 0, 5000, 10, 0)
        form.addRow("tophat radius (px)", self.tophat)
        form.addRow("opening radius (px)", self.popen)
        form.addRow("membrane halo excl. (px)", self.memex)
        form.addRow("fixed threshold", self.pthresh)
        self.btn_punc = QPushButton("Detect puncta"); self.btn_punc.clicked.connect(self._detect_puncta)
        form.addRow(self.btn_punc)
        self.punc_out = self._status(); form.addRow(self.punc_out)
        outer.addWidget(gb)

        # ---------- Measure / Save ----------
        gb, form = _group("Measure && save")
        self.results_edit = QLineEdit()
        self.results_edit.setPlaceholderText("results folder to update (e.g. results_all)…")
        rbrowse = QPushButton("Browse"); rbrowse.clicked.connect(self._browse_results)
        rrow = QHBoxLayout(); rrow.addWidget(self.results_edit); rrow.addWidget(rbrowse)
        form.addRow("Results folder", self._wrap(rrow))
        self.btn_table = QPushButton("Compute full metrics"); self.btn_table.clicked.connect(self._metrics)
        self.btn_save = QPushButton("Save / update results"); self.btn_save.clicked.connect(self._save)
        form.addRow(self.btn_table); form.addRow(self.btn_save)
        self.summary = self._status(); form.addRow(self.summary)
        outer.addWidget(gb)

        # ---------- Batch ----------
        gb, form = _group("Batch")
        self.folder_edit = QLineEdit(); self.folder_edit.setPlaceholderText("folder of ch1 TIFFs…")
        fbrowse = QPushButton("Browse"); fbrowse.clicked.connect(self._browse_folder)
        row = QHBoxLayout(); row.addWidget(self.folder_edit); row.addWidget(fbrowse)
        form.addRow("Folder", self._wrap(row))
        self.btn_batch = QPushButton("Batch process folder"); self.btn_batch.clicked.connect(self._batch)
        form.addRow(self.btn_batch)
        self.batch_prog = QProgressBar(); self.batch_prog.setFormat("%v / %m fields")
        self.batch_prog.setVisible(False)
        form.addRow(self.batch_prog)
        self.batch_out = self._status(); form.addRow(self.batch_out)
        outer.addWidget(gb)

        outer.addStretch(1)

    # ---------- small ui helpers ----------
    def _status(self):
        lab = QLabel(""); lab.setWordWrap(True)
        lab.setStyleSheet("color:#8ab4f8;")
        return lab

    def _wrap(self, layout):
        w = QWidget(); w.setLayout(layout); return w

    # ---------- engine helpers ----------
    def _ensure_maps(self):
        if self._bg is None:
            self._bg = compute_bg(self.image)

    def _seg_masks(self, img):
        d = self.diameter.value() or None
        m = segment_sam(img, diameter=d, flow_threshold=self.flow.value(),
                        cellprob_threshold=self.cellprob.value(), gpu=self.use_gpu.isChecked())
        m, _ = merge_fragments(m, min_area_frac=self.merge_frac.value())
        if self.rm_iso.isChecked():
            m, _ = remove_isolated(m)
        return m

    def _analyze(self, img, masks, use_cache=True):
        return analyze(img, masks, ring_px=self.ring.value(),
                       tophat_radius=self.tophat.value(), punc_open=self.popen.value(),
                       membrane_exclude_px=self.memex.value(),
                       punc_fixed_thresh=self.pthresh.value(),
                       bright_floor_k=self.bfloor.value(), bright_band_px=self.bband.value(),
                       bg=self._bg if use_cache else None)

    def _cur_masks(self):
        return np.asarray(self.viewer.layers["cells"].data).astype(np.int32)

    def _set_layer(self, name, data, kind, **kw):
        if name in self.viewer.layers:
            self.viewer.layers[name].data = data
        elif kind == "labels":
            self.viewer.add_labels(data, name=name, **kw)
        else:
            self.viewer.add_image(data, name=name, **kw)

    # ---------- callbacks ----------
    def _browse(self):
        f, _ = QFileDialog.getOpenFileName(self, "Filipin TIFF", "", "TIFF (*.tif *.tiff)")
        if f:
            self.path_edit.setText(f)

    def _browse_folder(self):
        d = QFileDialog.getExistingDirectory(self, "Batch folder")
        if d:
            self.folder_edit.setText(d)

    def _load(self):
        self.path = self.path_edit.text().strip()
        if not self.path or not os.path.exists(self.path):
            self.seg_out.setText("pick a valid TIFF"); return
        self.image = load_gray(self.path)
        self._bg = None; self.res = None
        # clear any layers from a previous image so overlays don't carry over
        for name in ("filipin", "cells", "membrane", "bright_membrane", "puncta"):
            if name in self.viewer.layers:
                self.viewer.layers.remove(name)
        for lab in (self.memb_out, self.punc_out, self.summary):
            lab.setText("")
        vmax = float(np.percentile(self.image, 99.5))
        self.viewer.add_image(self.image, name="filipin",
                              contrast_limits=[float(self.image.min()), vmax], colormap="gray")
        self.viewer.reset_view()
        self.seg_out.setText(f"loaded {os.path.basename(self.path)}  {self.image.shape}")

    def _segment(self):
        if self.image is None:
            self.seg_out.setText("load an image first"); return
        self.seg_out.setText("segmenting… (GPU ~15s, CPU minutes)"); self.repaint()
        masks = self._seg_masks(self.image)
        self._set_layer("cells", masks.astype(np.int32), "labels")
        self.viewer.layers["cells"].editable = True
        self.viewer.layers["cells"].opacity = 0.2
        self.seg_out.setText(f"{int(masks.max())} cells — edit 'cells' (paint/fill/erase), "
                             f"then Detect membrane / puncta.")

    def _detect_membrane(self):
        if self.image is None or "cells" not in self.viewer.layers:
            self.memb_out.setText("segment first"); return
        self.memb_out.setText("detecting membrane…"); self.repaint()
        self._ensure_maps()
        self.res = self._analyze(self.image, self._cur_masks())
        self._set_layer("membrane", self.res["membrane"], "image",
                        colormap="green", blending="additive", opacity=0.5)
        self._set_layer("bright_membrane", self.res["bright_membrane"], "image",
                        colormap="magenta", blending="additive", opacity=0.9)
        df = self.res["df"]
        self.memb_out.setText(
            f"ring membrane {df.membrane_mean.mean():.0f} (green) | bright membrane "
            f"cover {100*df.bright_membrane_cover.mean():.0f}% (magenta) | membrane/cyto "
            f"{df.membrane_to_cyto_ratio.mean():.2f} | membrane frac "
            f"{df.membrane_frac.mean():.2f}")

    def _detect_puncta(self):
        if self.image is None or "cells" not in self.viewer.layers:
            self.punc_out.setText("segment first"); return
        self.punc_out.setText("detecting puncta…"); self.repaint()
        self._ensure_maps()
        self.res = self._analyze(self.image, self._cur_masks())
        self._set_layer("puncta", self.res["pmask"], "image", colormap="cyan",
                        blending="additive", opacity=0.9)
        df = self.res["df"]
        self.punc_out.setText(
            f"puncta cover {100*df.puncta_cover_fixed.mean():.1f}% of cytoplasm | "
            f"puncta frac {df.puncta_frac.mean():.2f}  (cyan overlay)")

    def _metrics(self):
        if self.image is None or "cells" not in self.viewer.layers:
            self.summary.setText("segment first"); return
        self.summary.setText("computing…"); self.repaint()
        self._ensure_maps()
        self.res = self._analyze(self.image, self._cur_masks())
        df = self.res["df"]
        self.summary.setText(
            f"{len(df)} cells | membrane {df.membrane_mean.mean():.0f} "
            f"cyto {df.cytoplasm_mean.mean():.0f} ratio {df.membrane_to_cyto_ratio.mean():.2f} | "
            f"fracs M {df.membrane_frac.mean():.2f} C {df.cytoplasm_frac.mean():.2f} "
            f"P {df.puncta_frac.mean():.2f} | puncta {100*df.puncta_cover_fixed.mean():.1f}%")

    def _params(self):
        return dict(diameter=self.diameter.value(), flow=self.flow.value(),
                    cellprob=self.cellprob.value(), merge_frac=self.merge_frac.value(),
                    remove_isolated=self.rm_iso.isChecked(), use_gpu=self.use_gpu.isChecked(),
                    ring=self.ring.value(), bright_floor_k=self.bfloor.value(),
                    bright_band_px=self.bband.value(), tophat=self.tophat.value(),
                    opening=self.popen.value(), membrane_exclude=self.memex.value(),
                    punc_fixed_thresh=self.pthresh.value())

    def _browse_results(self):
        d = QFileDialog.getExistingDirectory(self, "Results folder")
        if d:
            self.results_edit.setText(d)

    def _results_dir(self):
        """Canonical results folder to write/update. Explicit field if set,
        else <batch folder>/filipin_results, else next to the image."""
        r = self.results_edit.text().strip()
        if r:
            return r
        fold = self.folder_edit.text().strip()
        if fold:
            return os.path.join(fold, "filipin_results")
        return os.path.join(os.path.dirname(self.path), "filipin_results")

    def _save(self):
        """Save/refine THIS field into the results folder and rebuild the aggregate,
        so a manual correction replaces this field's rows in all_cells.csv/summaries."""
        if self.res is None:
            self.summary.setText("compute metrics first"); return
        outdir = self._results_dir()
        stem = save_field(outdir, self.path, self._cur_masks(), self.res, img=self.image)
        n = rebuild_aggregate(outdir)
        with open(os.path.join(outdir, "params.json"), "w") as fh:
            json.dump(self._params(), fh, indent=2)
        self.summary.setText(f"updated '{stem}' in {os.path.basename(outdir)} "
                             f"→ rebuilt all_cells.csv ({n} fields)")

    def _batch(self):
        folder = self.folder_edit.text().strip()
        files = sorted(set(f for g in CH1_GLOBS
                           for f in glob.glob(os.path.join(folder, "**", g), recursive=True)))
        if not files:
            self.batch_out.setText(f"no {' / '.join(CH1_GLOBS)} under folder"); return
        outdir = self._results_dir()
        self.results_edit.setText(outdir)          # so later single Saves update here
        n_files = len(files)
        self.btn_batch.setEnabled(False)
        self.batch_prog.setRange(0, n_files); self.batch_prog.setValue(0)
        self.batch_prog.setVisible(True)
        failed = []
        try:
            for i, f in enumerate(files):
                name = os.path.basename(f)

                def show(msg):
                    self.batch_out.setText(f"[{i+1}/{n_files}] {name} — {msg}")
                    QApplication.processEvents()   # keep the window responsive

                try:
                    show("loading"); img = load_gray(f)
                    show("segmenting"); masks = self._seg_masks(img)
                    show("measuring"); res = self._analyze(img, masks, use_cache=False)
                    show("saving"); save_field(outdir, f, masks, res, img=img)
                except Exception as e:                       # one bad field mustn't abort the run
                    failed.append(name)
                    self.batch_out.setText(f"[{i+1}/{n_files}] {name} — FAILED: {e}")
                self.batch_prog.setValue(i + 1)
                QApplication.processEvents()

            n = rebuild_aggregate(outdir)
            with open(os.path.join(outdir, "params.json"), "w") as fh:
                json.dump(self._params(), fh, indent=2)
            msg = f"done: {n} fields → {os.path.basename(outdir)}/all_cells.csv"
            if failed:
                shown = ", ".join(failed[:3]) + ("…" if len(failed) > 3 else "")
                msg += f"  ({len(failed)} failed: {shown})"
            self.batch_out.setText(msg)
        finally:
            self.btn_batch.setEnabled(True)
            QApplication.processEvents()


def main():
    viewer = napari.Viewer(title="Filipin quantification")
    w = FilipinWidget(viewer)
    scroll = QScrollArea(); scroll.setWidgetResizable(True); scroll.setWidget(w)
    scroll.setMinimumWidth(300)
    viewer.window.add_dock_widget(scroll, area="right", name="Filipin quant")
    napari.run()


if __name__ == "__main__":
    main()
