"""Segment paired brightfield TIFFs, then score GFP inside reviewed cell masks."""

import argparse
import csv
import hashlib
from importlib.metadata import version
import json
import os
from pathlib import Path

import cv2
import numpy as np
from scipy import ndimage as ndi


def sample_name(path):
    # RLMI_like files end in a space-separated acquisition identifier.
    return path.stem.rsplit(" ", 1)[0]


def sample_key(sample):
    if sample.startswith("RLMI-") and sample[5:].isdigit():
        return (0, int(sample[5:]))
    return (1, sample.lower())


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def paired_images(root, requested):
    indexed = {}
    for channel in ("BRIGHTFIELD", "GREEN"):
        folder = root / channel
        files = sorted(folder.glob("*.tif"))
        if not files:
            raise ValueError(f"No TIFF files in {folder}")
        entries = {}
        for path in files:
            sample = sample_name(path)
            if sample in entries:
                raise ValueError(f"Duplicate sample {sample} in {folder}")
            entries[sample] = path
        indexed[channel] = entries
    if indexed["BRIGHTFIELD"].keys() != indexed["GREEN"].keys():
        raise ValueError("BRIGHTFIELD and GREEN sample labels do not match")
    if requested and requested not in indexed["GREEN"]:
        raise ValueError(f"Sample not found: {requested}")
    samples = [requested] if requested else sorted(indexed["GREEN"], key=sample_key)
    for sample in samples:
        if indexed["BRIGHTFIELD"][sample].stem != indexed["GREEN"][sample].stem:
            raise ValueError(f"Acquisition suffix differs between channels for {sample}")
    return [(sample, indexed["BRIGHTFIELD"][sample], indexed["GREEN"][sample]) for sample in samples]


def read_pair(brightfield_path, green_path, crop_bottom):
    brightfield = cv2.imread(str(brightfield_path), cv2.IMREAD_UNCHANGED)
    fluorescence = cv2.imread(str(green_path), cv2.IMREAD_UNCHANGED)
    if brightfield is None or fluorescence is None:
        raise ValueError(f"Unreadable image pair: {brightfield_path}, {green_path}")
    if brightfield.dtype != np.uint8 or brightfield.ndim != 2:
        raise ValueError(f"Expected 8-bit grayscale brightfield: {brightfield_path}")
    if fluorescence.dtype != np.uint8 or fluorescence.ndim != 3 or fluorescence.shape[2] < 3:
        raise ValueError(f"Expected 8-bit multichannel GREEN TIFF: {green_path}")
    if fluorescence.shape[:2] != brightfield.shape:
        raise ValueError(f"Channel dimensions differ: {brightfield_path}, {green_path}")
    if crop_bottom is not None and crop_bottom > brightfield.shape[0]:
        raise ValueError(f"Crop exceeds image height: {brightfield_path}")
    limit = crop_bottom or brightfield.shape[0]
    return brightfield[:limit], fluorescence[:limit, :, 1]


def write_csv(path, records):
    if not records:
        raise ValueError(f"No records for {path}")
    with path.open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(records[0]))
        writer.writeheader()
        writer.writerows(records)


def validate_masks(labels, shape, sample):
    if labels is None or labels.shape != shape or labels.dtype != np.uint16:
        raise ValueError(f"Invalid uint16 mask or shape for {sample}")
    count = int(labels.max())
    if count == 0 or not np.array_equal(np.unique(labels), np.arange(count + 1)):
        raise ValueError(f"Empty or non-contiguous mask labels for {sample}")
    return count


def overlay(brightfield, green, labels, positive=None):
    gray = cv2.normalize(brightfield, None, 25, 145, cv2.NORM_MINMAX)
    rgb = np.stack((gray, np.clip(gray.astype(np.float32) + green * .65, 0, 255), gray), axis=-1).astype(np.uint8)
    boundary = (labels > 0) & ((labels != np.roll(labels, 1, axis=0)) |
                               (labels != np.roll(labels, 1, axis=1)))
    if positive is None:
        rgb[boundary] = (0, 220, 255)
    else:
        selected = np.zeros(len(positive) + 1, dtype=bool)
        selected[1:] = positive
        rgb[boundary & ~selected[labels]] = (0, 220, 255)
        rgb[boundary & selected[labels]] = (255, 60, 200)
    return rgb


def save_qc(path, rgb, full):
    if not full:
        rgb = cv2.resize(rgb, (rgb.shape[1] // 2, rgb.shape[0] // 2), interpolation=cv2.INTER_AREA)
    if not cv2.imwrite(str(path), cv2.cvtColor(rgb, cv2.COLOR_RGB2BGR)):
        raise OSError(f"Could not save QC image: {path}")


def contact_sheet(output_dir, pairs, image_name, destination):
    tile_width, tile_height, header, columns = 432, 297, 30, 4
    rows = (len(pairs) + columns - 1) // columns
    sheet = np.full((rows * (tile_height + header), columns * tile_width, 3), 255, np.uint8)
    for index, (sample, _, _) in enumerate(pairs):
        image = cv2.imread(str(output_dir / image_name(sample)), cv2.IMREAD_COLOR)
        if image is None:
            raise ValueError(f"Missing QC overview for {sample}")
        x, y = index % columns * tile_width, index // columns * (tile_height + header)
        sheet[y + header:y + header + tile_height, x:x + tile_width] = cv2.resize(
            image, (tile_width, tile_height), interpolation=cv2.INTER_AREA)
        cv2.putText(sheet, sample, (x + 7, y + 21), cv2.FONT_HERSHEY_SIMPLEX,
                    .6, (0, 0, 0), 1, cv2.LINE_AA)
    if not cv2.imwrite(str(output_dir / destination), sheet, [cv2.IMWRITE_JPEG_QUALITY, 92]):
        raise OSError(f"Could not save contact sheet: {destination}")


def segment(args, pairs):
    os.environ.setdefault("CELLPOSE_LOCAL_MODELS_PATH", str(args.output_dir / ".cellpose_models"))
    from cellpose import models

    args.output_dir.mkdir(parents=True, exist_ok=True)
    manifest_path = args.output_dir / "mask_manifest.json"
    if manifest_path.exists() or any(args.output_dir.glob("mask_*.tif")):
        raise ValueError("Masks already exist; choose a fresh output directory")
    model = models.Cellpose(gpu=False, model_type=args.model)
    model_path = Path(model.cp.pretrained_model)
    manifest = {"model": args.model, "diameter": args.diameter, "invert": args.invert,
                "crop_bottom": args.crop_bottom, "cellpose_version": version("cellpose"),
                "model_sha256": sha256(model_path), "samples": {}}
    for sample, brightfield_path, green_path in pairs:
        brightfield, green = read_pair(brightfield_path, green_path, args.crop_bottom)
        labels, _, _, _ = model.eval(brightfield, channels=[0, 0], diameter=args.diameter, invert=args.invert)
        if int(labels.max()) > np.iinfo(np.uint16).max:
            raise ValueError(f"Too many labels for uint16 mask: {sample}")
        labels = labels.astype(np.uint16)
        count = validate_masks(labels, brightfield.shape, sample)
        mask_path = args.output_dir / f"mask_{sample}.tif"
        if not cv2.imwrite(str(mask_path), labels):
            raise OSError(f"Could not save mask: {mask_path}")
        save_qc(args.output_dir / f"qc_mask_{sample}_overview.png", overlay(brightfield, green, labels), False)
        if sample in args.full_qc:
            save_qc(args.output_dir / f"qc_mask_{sample}_full.png", overlay(brightfield, green, labels), True)
        manifest["samples"][sample] = {
            "brightfield_file": brightfield_path.name, "green_file": green_path.name,
            "brightfield_sha256": sha256(brightfield_path), "green_sha256": sha256(green_path),
            "mask_sha256": sha256(mask_path), "mask_count": count,
        }
        print(f"{sample}: {count} masks", flush=True)
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    contact_sheet(args.output_dir, pairs, lambda sample: f"qc_mask_{sample}_overview.png",
                  "qc_mask_all_samples_contact_sheet.jpg")


def verify_manifest(args, pairs, mask_dir):
    manifest_path = mask_dir / "mask_manifest.json"
    if not manifest_path.exists():
        if not args.accept_unverified_masks:
            raise ValueError("Mask manifest missing; use --accept-unverified-masks only for independently reviewed masks")
        return "external masks; source and model provenance not verified"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    if manifest["crop_bottom"] != args.crop_bottom:
        raise ValueError("Crop differs from mask manifest")
    for sample, brightfield_path, green_path in pairs:
        entry = manifest["samples"].get(sample)
        mask_path = mask_dir / f"mask_{sample}.tif"
        if (not entry or entry["brightfield_file"] != brightfield_path.name or
            entry["green_file"] != green_path.name or
            entry["brightfield_sha256"] != sha256(brightfield_path) or
            entry["green_sha256"] != sha256(green_path) or
            entry["mask_sha256"] != sha256(mask_path)):
            raise ValueError(f"Input/mask provenance mismatch for {sample}")
    return (f"{manifest['model']} {manifest['cellpose_version']}; "
            f"model SHA-256={manifest['model_sha256']}; diameter={manifest['diameter']}; invert={manifest['invert']}")


def score_cells(sample, labels, green, args):
    background = float(np.percentile(green, args.background_percentile))
    kernel = np.ones((2 * args.erode_pixels + 1,) * 2, dtype=np.uint8)
    scores, records = [], []
    for label, region in enumerate(ndi.find_objects(labels), start=1):
        if region is None:
            raise ValueError(f"Missing mask label {label}: {sample}")
        cell = (labels[region] == label).astype(np.uint8)
        inner = cv2.erode(cell, kernel, iterations=1).astype(bool) if args.erode_pixels else cell.astype(bool)
        if inner.sum() < 9:
            inner = cell.astype(bool)
        p = float(np.percentile(green[region][inner], args.cell_percentile))
        score = float(np.float32(p - background))
        scores.append(round(score, 3))
        records.append({"sample": sample, "mask_label": label, "mask_pixels": int(cell.sum()),
                        "interior_pixels": int(inner.sum()),
                        "gfp_interior_percentile_8bit": round(p, 3),
                        "image_background_percentile_8bit": round(background, 3),
                        "background_corrected_score_8bit": round(score, 3)})
    return np.array(scores), records


def borderline_plate(sample, labels, brightfield, green, scores, threshold, out):
    width = 112
    ys, xs = np.indices(labels.shape)
    weights = np.bincount(labels.ravel(), minlength=len(scores) + 1)
    cy = np.bincount(labels.ravel(), weights=ys.ravel(), minlength=len(scores) + 1) / np.maximum(weights, 1)
    cx = np.bincount(labels.ravel(), weights=xs.ravel(), minlength=len(scores) + 1) / np.maximum(weights, 1)
    candidates = np.flatnonzero((scores > threshold - 5) & (scores <= threshold)) + 1
    candidates = candidates[(cy[candidates] >= width // 2) & (cy[candidates] < labels.shape[0] - width // 2) &
                            (cx[candidates] >= width // 2) & (cx[candidates] < labels.shape[1] - width // 2)]
    if not len(candidates):
        return
    selected = candidates[np.linspace(0, len(candidates) - 1, min(12, len(candidates)), dtype=int)]
    plate = np.full(((width + 30) * 4, width * 6, 3), 255, np.uint8)
    for index, label in enumerate(selected):
        top, left = int(round(cy[label] - width / 2)), int(round(cx[label] - width / 2))
        section = (slice(top, top + width), slice(left, left + width))
        raw = green[section]
        gray = cv2.normalize(brightfield[section], None, 25, 145, cv2.NORM_MINMAX)
        annotated = np.stack((gray, np.clip(gray.astype(float) + raw * .65, 0, 255), gray), axis=-1).astype(np.uint8)
        instance = (labels[section] == label).astype(np.uint8)
        boundary = (instance > 0) & (cv2.erode(instance, np.ones((3, 3), np.uint8)) == 0)
        annotated[boundary] = (255, 210, 0)
        y, x = (index // 3) * (width + 30), (index % 3) * width * 2
        plate[y + 30:y + 30 + width, x:x + width] = np.stack((np.zeros_like(raw), raw, np.zeros_like(raw)), axis=-1)
        plate[y + 30:y + 30 + width, x + width:x + width * 2] = annotated
        cv2.putText(plate, f"{sample} #{label} S={scores[label - 1]:.1f}", (x + 3, y + 20),
                    cv2.FONT_HERSHEY_SIMPLEX, .5, (0, 0, 0), 1, cv2.LINE_AA)
    save_qc(out / f"qc_{sample}_borderline_T{threshold}.png", plate, True)


def classify(args, pairs):
    if not args.denominator_qc_reviewed:
        raise ValueError("Inspect full-resolution brightfield mask boundaries first, then pass --denominator-qc-reviewed")
    if (args.output_dir / "transfection_dual_index_report.csv").exists():
        raise ValueError("Classification already exists; choose a fresh output directory")
    mask_dir = args.mask_dir or args.output_dir
    provenance = verify_manifest(args, pairs, mask_dir)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    summary, all_cells, sensitivity, validation = [], [], [], []
    for sample, brightfield_path, green_path in pairs:
        brightfield, green = read_pair(brightfield_path, green_path, args.crop_bottom)
        labels = cv2.imread(str(mask_dir / f"mask_{sample}.tif"), cv2.IMREAD_UNCHANGED)
        count = validate_masks(labels, brightfield.shape, sample)
        scores, records = score_cells(sample, labels, green, args)
        if len(scores) != count:
            raise ValueError(f"Mask/score count mismatch for {sample}")
        positive = scores > args.threshold
        weights = np.clip((scores - args.threshold) / (255 - args.threshold), 0, 1)
        old = scores > 40
        rate = float(positive.mean() * 100)
        expression = float(weights.mean() * 100)
        mean_positive = float(weights[positive].mean() * 100) if positive.any() else 0.
        label = f"T{args.threshold}"
        summary.append({"sample": sample, "total_brightfield_masks": count,
                        f"gfp_positive_masks_{label}": int(positive.sum()),
                        f"estimated_cell_positive_percent_{label}": round(rate, 3),
                        f"expression_efficiency_percent_{label}": round(expression, 3),
                        f"mean_positive_relative_brightness_percent_{label}": round(mean_positive, 3),
                        "estimated_cell_positive_percent_T40": round(float(old.mean() * 100), 3),
                        "additional_positive_masks_vs_T40": int(positive.sum() - old.sum()),
                        "negative_control": "none", "false_positive_rate": "unknown",
                        "green_tiff": str(green_path.relative_to(args.input_dir)),
                        "brightfield_mask_tiff": str((mask_dir / f"mask_{sample}.tif").resolve())})
        for cell, yes, weight, legacy in zip(records, positive, weights, old):
            all_cells.append({"sample": sample, "mask_label": cell["mask_label"],
                              "background_corrected_score_8bit": cell["background_corrected_score_8bit"],
                              f"positive_at_{args.threshold}": int(yes), "positive_at_40": int(legacy),
                              f"relative_expression_weight_{label}": round(float(weight), 5)})
        write_csv(args.output_dir / f"cells_{sample}.csv", records)
        for threshold in args.sensitivity:
            sensitivity.append({"sample": sample, "threshold_above_image_background_8bit": threshold,
                                "positive_masks": int((scores > threshold).sum()), "total_masks": count,
                                "estimated_cell_positive_percent": round(float((scores > threshold).mean() * 100), 3)})
        validation.append({"sample": sample, "total_masks": count, "per_cell_records": len(records),
                           "expression_identity_residual_percent": round(abs(expression - rate * mean_positive / 100), 4),
                           "negative_control_false_positive_rate": "unknown"})
        rgb = overlay(brightfield, green, labels, positive)
        save_qc(args.output_dir / f"qc_{sample}_{label}_overview.png", rgb, False)
        if sample in args.full_qc:
            save_qc(args.output_dir / f"qc_{sample}_{label}_full.png", rgb, True)
            borderline_plate(sample, labels, brightfield, green, scores, args.threshold, args.output_dir)
        print(f"{sample}: {positive.sum()}/{count} positive at T={args.threshold}", flush=True)
    write_csv(args.output_dir / "transfection_dual_index_report.csv", summary)
    write_csv(args.output_dir / "per_cell_binary_classification.csv", all_cells)
    write_csv(args.output_dir / "threshold_sensitivity.csv", sensitivity)
    write_csv(args.output_dir / "method_validation.csv", validation)
    config = {"model_provenance": provenance, "threshold": args.threshold, "sensitivity": args.sensitivity,
              "crop_bottom": args.crop_bottom, "erode_pixels": args.erode_pixels,
              "cell_percentile": args.cell_percentile, "background_percentile": args.background_percentile,
              "denominator_qc": "operator-reviewed (not independently validated)",
              "negative_control": "none", "false_positive_rate": "unknown"}
    (args.output_dir / "analysis_parameters.json").write_text(json.dumps(config, indent=2) + "\n", encoding="utf-8")
    lines = [f"# Exploratory GFP cell-positive fraction (T={args.threshold}, uncalibrated)", "",
             f"{len(summary)} samples, {len(all_cells)} brightfield masks. Cellpose: {provenance}.",
             f"Score = eroded cell interior green p{args.cell_percentile:g} - image green p{args.background_percentile:g}; "
             f"8-bit TIFF, crop first {args.crop_bottom or 'all'} rows, erosion {args.erode_pixels} px; positive iff S > {args.threshold}.",
             "Denominator masks require manual QC; --denominator-qc-reviewed is a user assertion, not an automatic validation.",
             "Negative control: none; false-positive rate: unknown. Visually invisible background is an acquisition assumption, not a measurement of zero background.",
             "Exposure/gain/contrast changes can alter detection sensitivity and brightness; avoid strict cross-sample ranking without matched acquisition.",
             "Sensitivity cutoffs are not confidence intervals. Inspect full-resolution mask and borderline QC images.", "",
             f"| Sample | Positive % (T={args.threshold}) | Positive/total | Expression index % |", "|---|---:|---:|---:|"]
    for row in summary:
        lines.append(f"| {row['sample']} | {row[f'estimated_cell_positive_percent_T{args.threshold}']:.3f} | "
                     f"{row[f'gfp_positive_masks_T{args.threshold}']}/{row['total_brightfield_masks']} | "
                     f"{row[f'expression_efficiency_percent_T{args.threshold}']:.3f} |")
    lines.extend(["", "These are exploratory image-based indices, not validated absolute transfection efficiency.", ""])
    (args.output_dir / "transfection_dual_index_report.md").write_text("\n".join(lines), encoding="utf-8")
    contact_sheet(args.output_dir, pairs, lambda sample: f"qc_{sample}_T{args.threshold}_overview.png",
                  f"qc_all_samples_T{args.threshold}_contact_sheet.jpg")


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("stage", choices=("masks", "classify"))
    parser.add_argument("--input-dir", required=True, type=Path, help="Directory with BRIGHTFIELD/ and GREEN/ TIFFs")
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--sample", help="Optional single sample label for a trial run")
    parser.add_argument("--crop-bottom", type=int, default=None, help="Analyze only the first N rows")
    parser.add_argument("--model", default="cyto3", help="Cellpose 3 cytoplasm model for brightfield masks")
    parser.add_argument("--diameter", type=float, default=35.)
    parser.add_argument("--invert", action="store_true")
    parser.add_argument("--threshold", type=int, default=25)
    parser.add_argument("--sensitivity", type=int, nargs="+", default=(20, 25, 30, 40, 60))
    parser.add_argument("--cell-percentile", type=float, default=90.)
    parser.add_argument("--background-percentile", type=float, default=20.)
    parser.add_argument("--erode-pixels", type=int, default=2)
    parser.add_argument("--full-qc", nargs="*", default=(), help="Samples with full-size and numbered QC")
    parser.add_argument("--mask-dir", type=Path, help="Existing independently reviewed mask directory")
    parser.add_argument("--accept-unverified-masks", action="store_true", help="Allow legacy masks without a manifest")
    parser.add_argument("--denominator-qc-reviewed", action="store_true", help="Confirm manual mask-boundary review")
    args = parser.parse_args()
    if (args.crop_bottom is not None and args.crop_bottom <= 0 or args.diameter <= 0 or
        args.erode_pixels < 0 or args.threshold < 0 or args.threshold >= 255 or
        any(t < 0 or t > 255 for t in args.sensitivity) or
        any(p < 0 or p > 100 for p in (args.cell_percentile, args.background_percentile))):
        parser.error("Invalid crop, diameter, erosion, cutoff or percentile")
    return args


def main():
    args = parse_args()
    pairs = paired_images(args.input_dir, args.sample)
    if args.stage == "masks":
        segment(args, pairs)
    else:
        classify(args, pairs)


if __name__ == "__main__":
    main()
