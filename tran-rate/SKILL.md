---
name: tran-rate
description: "Estimate GFP or other reporter-positive cell fractions from paired microscopy images with manual cell-mask and threshold QC. Use for exploratory transfection comparisons with brightfield or nuclear images, including batches without a biological negative control; never claim validated absolute efficiency."
---

# Tran Rate

This is a simplified image-based estimation workflow, not a rigorous or
validated assay of absolute transfection efficiency. Use its outputs for rough,
exploratory comparison only when acquisition conditions are matched. Do not
present them as equivalent to flow cytometry or validated manual/nuclear cell
counting.

Preserve the original images and keep the calculation traceable. Do not use
fluorescence-positive object count as the default estimate; bright cells can
split into multiple objects and dim cells can disappear.

**Critical denominator QC gate:** Never accept inverted-brightfield peaks or
CellProfiler `TotalCells` simply because the pipeline ran or counts look
plausible. In a dense phase/brightfield RLMI_like batch, these peaks landed on
intercellular junctions; the resulting 24-sample rates were withdrawn. An
overlay with center markers is not sufficient: check actual cell interiors and
missed cells at full resolution before any percentage is published.

## Report two indices by default

For each fixed total-cell location (i), define:

- (S_i): local reporter signal after background correction.
- (T): a single-cell threshold calibrated from a matched negative control when
  available; otherwise an explicitly uncalibrated, manually reviewed cutoff.
- (M): maximum intensity used for normalization, such as 255 for 8-bit data.
- (C_i = 1) when (S_i > T), otherwise 0.
- (W_i = clamp((S_i - T) / (M - T), 0, 1)).

Report:

1. **Cell transfection rate (%)** = (mean(C_i) * 100). Each cell contributes
   only 0 or 1, regardless of how bright it is.
2. **Expression efficiency (%)** = (mean(W_i) * 100). This equals cell
   transfection rate multiplied by the mean relative brightness of positive
   cells, so it weights both prevalence and expression level.

Describe the first index as *not brightness-weighted*, not as physically
independent of fluorescence. Cells below the selected detection threshold
remain classified negative even when weak expression is visible on manual QC.

## Choose the total-cell denominator

- Prefer DAPI/Hoechst nuclei when present. Segment nuclei independently of the
  reporter channel.
- Use brightfield centers only when no nuclear channel exists. Treat the result
  as a relative estimate, especially in confluent cultures.
- Never segment the denominator from GFP when reporting cell transfection rate.
- For paired channels, derive cell masks solely from the original brightfield
  or nuclear TIFF and measure reporter signal separately within those masks
  in the original fluorescence TIFF. A merged image or JPG preview is for QC,
  not quantification. Confirm TIFF bit depth, channel mapping, pairing, and
  alignment; do not assume a display composite preserves raw values.
- For brightfield-only images, prefer an established segmentation tool such as
  Cellpose (try a cytoplasm model on raw and inverted brightfield), or a
  supervised ilastik pixel/object workflow with a small set of manually labeled
  examples. CellProfiler can pair images and measure *validated* masks. A
  pre-trained model is a starting point, not an automatic acceptance criterion.
- Inspect mask *boundaries* over raw brightfield and reporter images. Masks that
  partition cell-free gaps, merge adjacent cells, or miss dim/unusual cells
  invalidate the denominator. Correct representative masks in the tool and
  fine-tune/train if needed; otherwise stop at image-level relative fluorescence.
- At full resolution, inspect whether masks cover cell interiors rather than
  intercellular junctions. Compare low-, middle-, and high-signal fields;
  review dim cells, masks close to bright neighbors, and cell-free regions.

## Prepare and process the batch

1. Inventory original TIFF files and identify paired reporter and denominator
   images from the same field. Do not quantify JPG previews.
2. Standardize pairs without modifying originals:

   ```powershell
   python scripts/prepare_pairs.py --input-dir <images> --output-dir <paired> --alias o=0 --alias lipo2000=2000
   ```

3. Review `pair_manifest.csv`; stop for missing or duplicate channels.
4. Segment cells on the denominator channel with Cellpose or an annotated
   ilastik/CellProfiler workflow. The bundled
   `assets/transfection_area_and_count.cppipe` is an **unvalidated legacy
   starting point** for dense brightfield and must not be used to publish a
   denominator without manual QC. Its inverted-brightfield enhancement can
   detect cell junctions as centers. Treat `PositiveCells.csv` and GFP-object
   ratios as legacy diagnostics only.
5. Before classification, manually audit at least low-, medium-, and
   high-reporter fields, including both dense and sparse regions. Compare
   boundaries to cells in the *brightfield*, record clearly missed, duplicate,
   merged, and gap masks, and review at least 100 candidate cells per field
   plus a scan for missed cells. Keep annotated overlays and the audit counts.
   Reject the batch if a systematic error remains; do not repair a bad
   denominator by changing GFP thresholds.
6. Only after denominator QC passes, calculate both indices. For validated
   CellProfiler centers **with a biological negative control**, the bundled
   classifier is available:

   ```powershell
   python scripts/classify_cells_dual_index.py `
     --cells-csv <cellprofiler-output>/TotalCells.csv `
     --input-dir <paired> `
     --output-dir <dual-index-output> `
     --negative-control <sample>
   ```

   Add `--crop-bottom`, `--cell-radius`, `--negative-quantile`, or
   `--intensity-max` when the defaults do not match the batch. Use
   `--intensity-max 255` for 8-bit data and verify raw-value preservation
   before using a different bit depth. This bundled classifier **requires a
   biological negative control** and center coordinates; it is not a direct
   importer of Cellpose masks. With no control or with masks, use an appropriate
   mask-based measurement/classifier, inspect score locations on overlays,
   and label its threshold uncalibrated. Do not pass a transfected sample as
   `--negative-control` merely to make the bundled script run.

## Calibrate the classifier

- When available, use a biological negative control acquired with the same
  exposure, gain, magnification, and processing.
- The bundled script measures the green-channel 90th percentile within a fixed
  disk around each total-cell center, then subtracts the image background
  median.
- Default to the negative-control 99.5th percentile for (T), targeting about
  0.5% negative-control false positives. State the observed rate.
- Apply the same (T) to every image in the acquisition batch.
- Do not reuse (T), crop, radius, or intensity maximum across changed imaging
  settings without recalibration.
- With no biological negative control, a fixed heuristic background-corrected
  cutoff may be explored only after denominator QC. If the operator says
  background/autofluorescence was adjusted to be visually invisible, record
  that as an acquisition assumption, not proof that its pixels are zero or
  that the false-positive rate is zero. Never silently nominate a dim
  transfected image as the negative control.
- Choose an exploratory cutoff by comparing several thresholds and auditing
  numbered, full-resolution borderline cells side by side in the unadjusted
  reporter TIFF and brightfield overlay. Check both newly positive masks and
  visibly positive masks still below threshold; distinguish a sensitivity
  choice from proof of specificity. Report the sensitivity table, label rates
  uncalibrated, and state the false-positive rate is unknown. The threshold
  range is not a confidence interval.
- A background-subtracted cutoff in 8-bit gray levels is a fixed number *in
  the stored image*, not an absolute fluorescence quantity. Exposure, gain,
  contrast, bit depth, or processing can alter the score even after subtracting
  each image's background. Record the exact score formula and strict `S > T`
  rule. Do not carry a numeric cutoff across changed acquisition settings
  without recalibration and new manual QC; per-image background subtraction
  alone does not normalize sensitivity across samples.
- RLMI_like case example, not a default: interior green-channel p90 minus
  image-wide green-channel p20 in 8-bit original TIFFs; a 2-pixel mask erosion
  and `T=25` were a user-selected exploratory tradeoff after reviewing `T=20`
  borderline cells. Many cells scoring 20-25 appeared positive, so `T=25`
  knowingly leaves some false negatives. This number must not be hard-coded
  into other datasets or presented as an optimal calibrated threshold.

## Prevent neighboring-cell contamination

Test several sampling radii appropriate to cell size, for example 4, 6, 8, 10,
and 12 pixels. If the positive rate rises monotonically with radius in
confluent fields, reporter signal from neighboring cells is entering the
measurement. Select the smallest radius that still samples the intended
intracellular region, and record the tested range and chosen value.

## Validate before accepting results

- Confirm one reporter and one denominator image per sample.
- Require manual mask/center-location QC to pass *before* calculating rates.
  Archive overlays at full resolution and a record of audited objects/errors.
- Confirm every CellProfiler `ModuleError_*` value is zero.
- Confirm unique `ImageNumber,ObjectNumber` pairs and no missing cell scores.
- Confirm sample total-cell counts equal the counts in `TotalCells.csv`.
- Inspect composite overlays for a negative control, a median sample, and a
  high-positive sample (or low-, median-, and high-signal samples if no
  negative exists). Green is raw reporter signal. For a mask workflow, draw
  actual cell boundaries rather than only center dots; count obvious missed
  cells and false masks on dark intercellular areas.
- Show raw reporter crops next to full-resolution boundaries for borderline
  scores; overlays with enhanced display contrast alone can conceal faint
  positives or neighboring-cell spillover.
- When a negative control exists, review its 99%, 99.5%, and 99.9%
  threshold sensitivity table.
- Verify:

  ```text
  expression efficiency
    = cell transfection rate * mean positive-cell relative brightness / 100
  ```

- Flag saturated images. Expression efficiency is comparable only under the
  same acquisition settings. If exposure/gain/contrast may vary among samples,
  even the binary positive fractions can have different detection limits;
  avoid a strict cross-sample ranking without matched acquisition or controls.

## Required outputs after denominator QC passes

Keep these files in a fresh output folder:

- `transfection_dual_index_report.csv`: concise final two-index table.
- `transfection_dual_index_report.md`: definitions, results, and limitations.
- `per_cell_binary_classification.csv`: auditable per-cell scores and labels.
- `threshold_sensitivity.csv`: cutoff sensitivity.
- `method_validation.csv`: false-positive and correlation checks.
- `qc_*.png`: negative, median, and high-positive overlays.

State the negative control, (T/M), cell radius, crop, denominator channel,
sample count, false-positive rate, and acquisition-comparability limit.
If no control exists, explicitly state "none" and "unknown" for control and
false-positive rate. If denominator QC fails, output a failure/QC report only,
not `transfection_dual_index_report.csv` or a ranking.

When brightfield supplies the denominator, include:

> Estimated from brightfield cell centers without nuclear stain; suitable for
> relative comparison, not a validated absolute transfection efficiency.

## Resources

- `scripts/prepare_pairs.py`: pair and normalize TIFF filenames safely.
- `assets/transfection_area_and_count.cppipe`: obtain brightfield total-cell
  centers and legacy diagnostics with CellProfiler.
- `scripts/classify_cells_dual_index.py`: generate the dual-index report and
  cell-level audit outputs.
- `scripts/summarize_cellprofiler.py`: summarize legacy area/object metrics
  only when comparison with earlier analyses is required.
