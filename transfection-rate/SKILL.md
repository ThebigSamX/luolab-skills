---
name: transfection-rate
description: "Analyze GFP or other reporter fluorescence microscopy batches and report two complementary indices by default: binary cell transfection rate (green versus not green, without brightness weighting) and brightness-weighted expression efficiency. Use for transfection-rate, transfection-efficiency, reporter-expression, fluorescence-positive cell counting, paired fluorescence plus brightfield images, DAPI/nuclear denominators, CellProfiler outputs, negative-control threshold calibration, or Chinese requests mentioning 转染率, 转染效率, 阳性细胞率, or 表达效率."
---

# Transfection Rate

Preserve the original images and produce an auditable dual-index report. Do not
use fluorescence-positive object count as the default cell transfection rate;
bright cells can split into multiple objects and dim cells can disappear.

## Report two indices by default

For each fixed total-cell location (i), define:

- (S_i): local reporter signal after background correction.
- (T): one shared single-cell threshold derived from negative controls.
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
independent of fluorescence. Reporter expression below the negative-control
detection distribution remains indistinguishable from a negative cell.

## Choose the total-cell denominator

- Prefer DAPI/Hoechst nuclei when present. Segment nuclei independently of the
  reporter channel.
- Use brightfield centers only when no nuclear channel exists. Treat the result
  as a relative estimate, especially in confluent cultures.
- Never segment the denominator from GFP when reporting cell transfection rate.

## Prepare and process the batch

1. Inventory original TIFF files and identify paired reporter and denominator
   images from the same field. Do not quantify JPG previews.
2. Standardize pairs without modifying originals:

   ```powershell
   python scripts/prepare_pairs.py --input-dir <images> --output-dir <paired> --alias o=0 --alias lipo2000=2000
   ```

3. Review `pair_manifest.csv`; stop for missing or duplicate channels.
4. Generate fixed total-cell locations. For the bundled brightfield workflow,
   import `assets/transfection_area_and_count.cppipe` into CellProfiler 4.2.x,
   run the standardized pairs, and use `TotalCells.csv`. Treat its
   `PositiveCells.csv` and GFP-object ratio as legacy diagnostics only.
5. Calculate both indices:

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
   before using a different bit depth.

## Calibrate the classifier

- Use a biological negative control acquired with the same exposure, gain,
  magnification, and processing.
- The bundled script measures the green-channel 90th percentile within a fixed
  disk around each total-cell center, then subtracts the image background
  median.
- Default to the negative-control 99.5th percentile for (T), targeting about
  0.5% negative-control false positives. State the observed rate.
- Apply the same (T) to every image in the acquisition batch.
- Do not reuse (T), crop, radius, or intensity maximum across changed imaging
  settings without recalibration.

## Prevent neighboring-cell contamination

Test several sampling radii appropriate to cell size, for example 4, 6, 8, 10,
and 12 pixels. If the positive rate rises monotonically with radius in
confluent fields, reporter signal from neighboring cells is entering the
measurement. Select the smallest radius that still samples the intended
intracellular region, and record the tested range and chosen value.

## Validate before accepting results

- Confirm one reporter and one denominator image per sample.
- Confirm every CellProfiler `ModuleError_*` value is zero.
- Confirm unique `ImageNumber,ObjectNumber` pairs and no missing cell scores.
- Confirm sample total-cell counts equal the counts in `TotalCells.csv`.
- Inspect composite overlays for a negative control, a median sample, and a
  high-positive sample. Green is raw reporter signal, cyan marks all total-cell
  centers, and magenta rings mark classified positives.
- Review the 99%, 99.5%, and 99.9% negative-control threshold sensitivity table.
- Verify:

  ```text
  expression efficiency
    = cell transfection rate * mean positive-cell relative brightness / 100
  ```

- Flag saturated images. Expression efficiency is comparable only under the
  same acquisition settings.

## Required outputs

Keep these files in a fresh output folder:

- `transfection_dual_index_report.csv`: concise final two-index table.
- `transfection_dual_index_report.md`: definitions, results, and limitations.
- `per_cell_binary_classification.csv`: auditable per-cell scores and labels.
- `threshold_sensitivity.csv`: cutoff sensitivity.
- `method_validation.csv`: false-positive and correlation checks.
- `qc_*.png`: negative, median, and high-positive overlays.

State the negative control, (T/M), cell radius, crop, denominator channel,
sample count, false-positive rate, and acquisition-comparability limit.

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
