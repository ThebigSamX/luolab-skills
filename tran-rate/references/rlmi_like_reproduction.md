# RLMI_like reproduction (no negative control)

This recipe reproduces the exploratory T=25 image-based indices from the 24
original paired 8-bit TIFF fields. Images and reviewed masks are **not** part
of the skill. Exact numerical reproduction needs those same source images,
matching Cellpose model weights and software, and a fresh manual review of
cell boundaries. Do not represent a regenerated denominator as independently
validated by the original user's QC.

## Inputs and environment

- Arrange original images under `<batch>/BRIGHTFIELD/` and `<batch>/GREEN/`.
  The files for a sample must have the same stem; the last space-separated
  part is an acquisition suffix, e.g. `RLMI-9 923.tif`. MERGE/JPG images
  are never used for quantification.
- Use Python 3.10 and Cellpose 3.1.1.3 `cyto3` on CPU. The verified RLMI_like
  environment used NumPy 2.0.2, SciPy 1.15.3, OpenCV headless 5.0.0.93,
  and PyTorch 2.14.0. Install matching versions when exact reproduction is
  needed; Cellpose may fetch model weights on first use. The cached cyto3
  model in the verified environment has SHA-256
  `2DC3087A8ABD7DA46D1AB0DDD5824639933CC3FF63B382AF3FA1939A392DB93C`.
- In an isolated Python 3.10 environment, install the verified packages with
  `python -m pip install "cellpose==3.1.1.3" "numpy==2.0.2" "scipy==1.15.3" "opencv-python-headless==5.0.0.93" "torch==2.14.0"`.
  Match the Cellpose model SHA-256 in the generated manifest; if this wheel
  combination is unavailable on your platform, record the versions actually
  used and compare resulting mask hashes/counts rather than claiming exact
  binary reproduction.
- The script rejects non-8-bit or misaligned image pairs, missing/duplicate
  sample labels, non-contiguous masks, and changed source TIFF hashes after
  mask creation. Keep output outside the source TIFF directories.

## Two stages

Run commands from the installed `tran-rate` skill directory (where `scripts/`
is located), with the dependencies installed. Replace `<batch>` and `<result>`
with absolute paths. Use a fresh output directory for each run:

```powershell
python scripts/cellpose_rate.py masks --input-dir <batch> --output-dir <result> `
  --crop-bottom 1780 --model cyto3 --diameter 35 `
  --full-qc RLMI-9 RLMI-17 RLMI-24 RLMI-14
```

Inspect `qc_mask_all_samples_contact_sheet.jpg`, `qc_mask_*_overview.png`,
and full-resolution `qc_mask_*_full.png` on
both cell interiors and cell-free junctions. Compare representative low-,
medium-, high-signal and dense/sparse regions. Record missed, merged, and
spurious cells. Stop if denominator QC is not acceptable; do not run the
classifier to repair a bad mask with a different GFP cutoff.

```powershell
python scripts/cellpose_rate.py classify --input-dir <batch> --output-dir <result> `
  --crop-bottom 1780 --threshold 25 --erode-pixels 2 `
  --cell-percentile 90 --background-percentile 20 `
  --denominator-qc-reviewed --full-qc RLMI-9 RLMI-17 RLMI-24 RLMI-14
```

The classifier uses only GREEN's stored G channel inside Cellpose masks.
Within each eroded mask it takes p90, subtracts the image-wide G p20, and
classifies strictly `S > 25`. It writes mask provenance (`mask_manifest.json`),
parameters, per-cell scores, a 20/25/30/40/60 sensitivity table, two-index
summary, validation table, all-sample QC contact sheet, full-size QC for
requested samples, and numbered
20-25 score review plates. The raw GREEN crop is on each plate's left;
yellow outlines mark target cells on the brightfield overlay at right.

In the verified analysis, the 24 samples had 65,041 masks. Example T=25
values: RLMI-9 = 612/2853 (21.451%), RLMI-17 = 498/2986 (16.678%),
RLMI-24 = 339/3023 (11.214%). Compare the full output table to the
original analysis when its CSV is available. The verified mask TIFF hashes
were RLMI-9 = `E7571DF37E82479F32FDA3C68455A5F47D2A038BF805BEECB09C9F8699C94E0B`
and RLMI-24 = `C7659E124F25BDAF3AFF7778AA7930D6B27197C967C4CF09B7E516371F7577BC`.
`--sample RLMI-24` may be used
for a trial in a separate output directory before processing the full batch.

`T=25` is an uncalibrated 8-bit cutoff selected for this dataset after manual
QC, not a universal fluorescence standard. Some 20-25 score cells were also
visibly positive. The user's assumption that background was visually hidden
at acquisition does not establish zero false positives. Different exposure,
gain, or contrast across samples can change detection sensitivity and makes
strict cross-sample ranking inappropriate.

The optional `--mask-dir` permits already reviewed external masks for a
classification-only comparison. If those masks have no manifest, the script
requires `--accept-unverified-masks` and reports their provenance as
unverified; that route does **not** demonstrate raw-to-mask reproduction.
