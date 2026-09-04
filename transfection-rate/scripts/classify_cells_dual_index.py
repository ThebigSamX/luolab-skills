"""Classify fixed brightfield-detected cells as GFP positive or negative.

This analysis reuses CellProfiler's TotalCells centers, measures GFP in the
same fixed-radius neighborhood around every center, and calibrates one binary
cutoff from a negative-control sample. A positive cell contributes exactly one
count regardless of how bright or large its GFP signal is.
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path

import numpy as np
import pandas as pd
from PIL import Image, ImageDraw, ImageOps


SAMPLE_PATTERN = re.compile(r"^\d+_(?P<sample>.+)_g\.tif$", re.IGNORECASE)


def read_cellprofiler_csv(path: Path) -> pd.DataFrame:
    """Read CellProfiler CSVs that may contain legacy Windows paths."""
    last_error: UnicodeDecodeError | None = None
    for encoding in ("utf-8-sig", "gb18030"):
        try:
            return pd.read_csv(path, encoding=encoding)
        except UnicodeDecodeError as error:
            last_error = error
    assert last_error is not None
    raise last_error


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Control-calibrated binary GFP classification per fixed cell."
    )
    parser.add_argument(
        "--cells-csv",
        type=Path,
        default=Path("cellprofiler_results_count/TotalCells.csv"),
    )
    parser.add_argument(
        "--input-dir", type=Path, default=Path("cell_count_input")
    )
    parser.add_argument(
        "--old-summary",
        type=Path,
        default=Path(
            "cellprofiler_results_count/"
            "transfection_area_and_cell_count_summary.csv"
        ),
    )
    parser.add_argument(
        "--output-dir", type=Path, default=Path("cell_count_binary_results")
    )
    parser.add_argument("--negative-control", default="2000")
    parser.add_argument("--crop-bottom", type=int, default=1780)
    parser.add_argument("--cell-radius", type=int, default=4)
    parser.add_argument("--negative-quantile", type=float, default=0.995)
    parser.add_argument(
        "--intensity-max",
        type=float,
        default=255.0,
        help="Maximum green-channel value for expression-efficiency normalization.",
    )
    return parser.parse_args()


def sample_from_filename(filename: str) -> str:
    match = SAMPLE_PATTERN.match(filename)
    if not match:
        raise ValueError(f"Cannot extract sample from GFP filename: {filename}")
    return match.group("sample")


def green_channel(path: Path, crop_bottom: int) -> np.ndarray:
    with Image.open(path) as image:
        array = np.asarray(image)
    if array.ndim == 2:
        green = array
    elif array.ndim == 3 and array.shape[2] >= 2:
        green = array[:, :, 1]
    else:
        raise ValueError(f"Unsupported image shape {array.shape} for {path}")
    return green[:crop_bottom].astype(np.float32, copy=False)


def disk_offsets(radius: int) -> tuple[np.ndarray, np.ndarray]:
    yy, xx = np.mgrid[-radius : radius + 1, -radius : radius + 1]
    mask = xx * xx + yy * yy <= radius * radius
    return yy[mask], xx[mask]


def measure_cells(
    cells: pd.DataFrame,
    image: np.ndarray,
    radius: int,
) -> pd.DataFrame:
    offset_y, offset_x = disk_offsets(radius)
    height, width = image.shape
    background = float(np.median(image))
    records: list[dict[str, float | int]] = []

    for row in cells.itertuples(index=False):
        center_x = int(round(float(row.Location_Center_X)))
        center_y = int(round(float(row.Location_Center_Y)))
        xs = center_x + offset_x
        ys = center_y + offset_y
        valid = (xs >= 0) & (xs < width) & (ys >= 0) & (ys < height)
        pixels = image[ys[valid], xs[valid]]
        if pixels.size == 0:
            raise ValueError(
                f"Cell {row.ObjectNumber} in image {row.ImageNumber} "
                "has no pixels inside the analysis region"
            )
        p90 = float(np.percentile(pixels, 90))
        records.append(
            {
                "ImageNumber": int(row.ImageNumber),
                "ObjectNumber": int(row.ObjectNumber),
                "center_x": float(row.Location_Center_X),
                "center_y": float(row.Location_Center_Y),
                "gfp_local_mean_0_to_255": float(np.mean(pixels)),
                "gfp_local_median_0_to_255": float(np.median(pixels)),
                "gfp_local_p90_0_to_255": p90,
                "gfp_local_max_0_to_255": float(np.max(pixels)),
                "image_background_median_0_to_255": background,
                "gfp_background_corrected_p90_0_to_255": p90 - background,
            }
        )
    return pd.DataFrame.from_records(records)


def render_qc(
    brightfield_path: Path,
    gfp_path: Path,
    cells: pd.DataFrame,
    destination: Path,
    crop_bottom: int,
) -> None:
    with Image.open(brightfield_path) as source:
        brightfield = source.convert("L").crop(
            (0, 0, source.width, crop_bottom)
        )
    brightfield_array = np.asarray(
        ImageOps.autocontrast(brightfield), dtype=np.float32
    )
    gfp_array = green_channel(gfp_path, crop_bottom)
    gray = brightfield_array * 0.42
    composite = np.stack(
        [
            gray,
            np.clip(gray + gfp_array * 0.85, 0, 255),
            gray,
        ],
        axis=2,
    ).astype(np.uint8)
    image = Image.fromarray(composite, mode="RGB")
    scale = 0.5
    image = image.resize(
        (round(image.width * scale), round(image.height * scale)),
        Image.Resampling.LANCZOS,
    )
    draw = ImageDraw.Draw(image)
    for row in cells.itertuples(index=False):
        x = float(row.center_x) * scale
        y = float(row.center_y) * scale
        draw.ellipse((x - 1, y - 1, x + 1, y + 1), fill=(0, 220, 255))
    for row in cells[cells["binary_gfp_positive"]].itertuples(index=False):
        x = float(row.center_x) * scale
        y = float(row.center_y) * scale
        draw.ellipse(
            (x - 3, y - 3, x + 3, y + 3),
            outline=(255, 80, 255),
            width=2,
        )
    image.save(destination)


def main() -> None:
    args = parse_args()
    if not 0 < args.negative_quantile < 1:
        raise ValueError("--negative-quantile must be between 0 and 1")
    if args.cell_radius < 1:
        raise ValueError("--cell-radius must be positive")
    if args.intensity_max <= 0:
        raise ValueError("--intensity-max must be positive")

    cells = read_cellprofiler_csv(args.cells_csv)
    required = {
        "ImageNumber",
        "ObjectNumber",
        "FileName_Brightfield",
        "FileName_GFPColor",
        "Location_Center_X",
        "Location_Center_Y",
    }
    missing = required.difference(cells.columns)
    if missing:
        raise ValueError(f"Missing TotalCells columns: {sorted(missing)}")

    args.output_dir.mkdir(parents=True, exist_ok=True)
    measured_groups: list[pd.DataFrame] = []
    image_metadata: list[dict[str, str | int]] = []

    grouped = cells.groupby("ImageNumber", sort=True)
    for image_number, group in grouped:
        gfp_names = group["FileName_GFPColor"].drop_duplicates().tolist()
        brightfield_names = group["FileName_Brightfield"].drop_duplicates().tolist()
        if len(gfp_names) != 1 or len(brightfield_names) != 1:
            raise ValueError(f"Image {image_number} does not have one paired filename")
        gfp_name = str(gfp_names[0])
        brightfield_name = str(brightfield_names[0])
        sample = sample_from_filename(gfp_name)
        gfp_path = args.input_dir / gfp_name
        brightfield_path = args.input_dir / brightfield_name
        if not gfp_path.is_file() or not brightfield_path.is_file():
            raise FileNotFoundError(f"Missing pair: {gfp_path}, {brightfield_path}")

        measured = measure_cells(
            group,
            green_channel(gfp_path, args.crop_bottom),
            args.cell_radius,
        )
        measured.insert(2, "sample", sample)
        measured.insert(3, "gfp_tif", gfp_name)
        measured.insert(4, "brightfield_tif", brightfield_name)
        measured_groups.append(measured)
        image_metadata.append(
            {
                "ImageNumber": int(image_number),
                "sample": sample,
                "gfp_tif": gfp_name,
                "brightfield_tif": brightfield_name,
            }
        )

    per_cell = pd.concat(measured_groups, ignore_index=True)
    negative = per_cell[per_cell["sample"].astype(str) == str(args.negative_control)]
    if negative.empty:
        raise ValueError(f"Negative control {args.negative_control!r} was not found")
    score_column = "gfp_background_corrected_p90_0_to_255"
    threshold = float(
        np.quantile(
            negative[score_column].to_numpy(),
            args.negative_quantile,
            method="higher",
        )
    )
    if threshold >= args.intensity_max:
        raise ValueError(
            "The negative-control threshold must be below --intensity-max"
        )
    per_cell["classification_threshold_0_to_255"] = threshold
    per_cell["binary_gfp_positive"] = per_cell[score_column] > threshold
    per_cell["gfp_threshold_excess_0_to_255"] = np.clip(
        per_cell[score_column] - threshold, 0.0, None
    )
    per_cell["relative_expression_weight_0_to_1"] = np.clip(
        per_cell["gfp_threshold_excess_0_to_255"]
        / (args.intensity_max - threshold),
        0.0,
        1.0,
    )

    summary = (
        per_cell.groupby(
            ["ImageNumber", "sample", "gfp_tif", "brightfield_tif"],
            sort=True,
            as_index=False,
        )
        .agg(
            total_brightfield_cells=("ObjectNumber", "count"),
            binary_gfp_positive_cells=("binary_gfp_positive", "sum"),
            summed_relative_expression_weight=(
                "relative_expression_weight_0_to_1",
                "sum",
            ),
            mean_cell_gfp_score_0_to_255=(score_column, "mean"),
            median_cell_gfp_score_0_to_255=(score_column, "median"),
        )
    )
    summary["binary_cell_count_transfection_percent"] = (
        summary["binary_gfp_positive_cells"]
        / summary["total_brightfield_cells"]
        * 100.0
    )
    summary["expression_brightness_weighted_efficiency_percent"] = (
        summary["summed_relative_expression_weight"]
        / summary["total_brightfield_cells"]
        * 100.0
    )
    summary["mean_relative_brightness_of_positive_cells_percent"] = np.where(
        summary["binary_gfp_positive_cells"] > 0,
        summary["summed_relative_expression_weight"]
        / summary["binary_gfp_positive_cells"]
        * 100.0,
        0.0,
    )
    summary["rank_by_binary_cell_count_ratio"] = (
        summary["binary_cell_count_transfection_percent"]
        .rank(method="min", ascending=False)
        .astype(int)
    )
    summary["rank_by_expression_efficiency"] = (
        summary["expression_brightness_weighted_efficiency_percent"]
        .rank(method="min", ascending=False)
        .astype(int)
    )
    summary["negative_control"] = str(args.negative_control)
    summary["classifier"] = "cell_local_p90_minus_image_background_median"
    summary["classification_threshold_0_to_255"] = threshold
    summary["negative_control_quantile"] = args.negative_quantile
    summary["intensity_max_for_normalization"] = args.intensity_max
    summary["cell_sampling_radius_pixels"] = args.cell_radius
    summary["crop_rows"] = f"0:{args.crop_bottom}"
    summary["interpretation"] = (
        "Binary count from fixed brightfield cell centers; one count per "
        "control-thresholded cell; relative estimate without nuclear stain"
    )

    if args.old_summary.is_file():
        old = read_cellprofiler_csv(args.old_summary).rename(
            columns={
                "estimated_cell_count_ratio_percent": "old_gfp_object_count_percent",
                "gfp_positive_area_percent": "old_gfp_positive_area_percent",
            }
        )
        old["sample"] = old["sample"].astype(str)
        summary["sample"] = summary["sample"].astype(str)
        summary = summary.merge(
            old[
                [
                    "sample",
                    "old_gfp_object_count_percent",
                    "old_gfp_positive_area_percent",
                ]
            ],
            on="sample",
            how="left",
            validate="one_to_one",
        )

    column_order = [
        "sample",
        "rank_by_binary_cell_count_ratio",
        "rank_by_expression_efficiency",
        "total_brightfield_cells",
        "binary_gfp_positive_cells",
        "binary_cell_count_transfection_percent",
        "expression_brightness_weighted_efficiency_percent",
        "mean_relative_brightness_of_positive_cells_percent",
        "old_gfp_object_count_percent",
        "old_gfp_positive_area_percent",
        "mean_cell_gfp_score_0_to_255",
        "median_cell_gfp_score_0_to_255",
        "classification_threshold_0_to_255",
        "negative_control",
        "negative_control_quantile",
        "intensity_max_for_normalization",
        "classifier",
        "cell_sampling_radius_pixels",
        "crop_rows",
        "gfp_tif",
        "brightfield_tif",
        "interpretation",
        "ImageNumber",
    ]
    column_order = [column for column in column_order if column in summary.columns]
    summary = summary[column_order].sort_values("ImageNumber")

    summary = summary.round(
        {
            "binary_cell_count_transfection_percent": 3,
            "expression_brightness_weighted_efficiency_percent": 3,
            "mean_relative_brightness_of_positive_cells_percent": 3,
            "summed_relative_expression_weight": 3,
            "old_gfp_object_count_percent": 3,
            "old_gfp_positive_area_percent": 3,
            "mean_cell_gfp_score_0_to_255": 3,
            "median_cell_gfp_score_0_to_255": 3,
            "classification_threshold_0_to_255": 3,
        }
    )
    per_cell_export = per_cell.copy()
    per_cell_export["binary_gfp_positive"] = per_cell_export[
        "binary_gfp_positive"
    ].astype(int)
    per_cell_export = per_cell_export.round(3)
    per_cell_export.to_csv(
        args.output_dir / "per_cell_binary_classification.csv", index=False
    )
    summary.to_csv(args.output_dir / "binary_cell_count_summary.csv", index=False)

    report_columns = [
        "sample",
        "binary_cell_count_transfection_percent",
        "expression_brightness_weighted_efficiency_percent",
        "mean_relative_brightness_of_positive_cells_percent",
        "binary_gfp_positive_cells",
        "total_brightfield_cells",
        "rank_by_binary_cell_count_ratio",
        "rank_by_expression_efficiency",
    ]
    dual_report = summary[report_columns].rename(
        columns={
            "binary_cell_count_transfection_percent": "cell_transfection_rate_percent",
            "expression_brightness_weighted_efficiency_percent": "expression_efficiency_percent",
            "mean_relative_brightness_of_positive_cells_percent": (
                "mean_positive_cell_relative_brightness_percent"
            ),
            "binary_gfp_positive_cells": "gfp_positive_cells",
            "total_brightfield_cells": "total_cells",
            "rank_by_binary_cell_count_ratio": "cell_transfection_rate_rank",
            "rank_by_expression_efficiency": "expression_efficiency_rank",
        }
    )
    dual_report.to_csv(
        args.output_dir / "transfection_dual_index_report.csv", index=False
    )
    report_lines = [
        "# 转染分析双指标报告",
        "",
        "## 指标定义",
        "",
        "1. **细胞转染率（%）**：GFP 阳性细胞数 / 明场识别总细胞数 × 100。每个细胞只按绿或不绿计 1/0，不按表达亮度加权。",
        "2. **表达效率（%）**：所有细胞的阈值以上 GFP 相对亮度之和 / 总细胞数 × 100。等价于细胞转染率 × 阳性细胞平均相对亮度。",
        "",
        f"阴性对照：`{args.negative_control}`；单细胞阳性阈值：`{threshold:.3f}/{args.intensity_max:g}`；细胞中心取样半径：`{args.cell_radius} px`；分析裁剪：`0:{args.crop_bottom}`。",
        "",
        "## 结果",
        "",
        "| 样本 | 细胞转染率（%） | 表达效率（%） | 阳性细胞平均相对亮度（%） | 阳性/总细胞 | 转染率排名 | 表达效率排名 |",
        "|---|---:|---:|---:|---:|---:|---:|",
    ]
    for row in dual_report.itertuples(index=False):
        report_lines.append(
            f"| {row.sample} | {row.cell_transfection_rate_percent:.3f} | "
            f"{row.expression_efficiency_percent:.3f} | "
            f"{row.mean_positive_cell_relative_brightness_percent:.3f} | "
            f"{row.gfp_positive_cells}/{row.total_cells} | "
            f"{row.cell_transfection_rate_rank} | {row.expression_efficiency_rank} |"
        )
    report_lines.extend(
        [
            "",
            "## 解释限制",
            "",
            "细胞转染率不按 GFP 亮度加权，但低于阴性背景与相机检测限的表达仍无法识别。表达效率受曝光、增益、饱和和位深影响，只能在相同成像设置的本批次内直接比较。总细胞数来自无核染条件下的明场细胞中心估算，适合相对比较，不代表经核染验证的绝对转染效率。",
            "",
        ]
    )
    (args.output_dir / "transfection_dual_index_report.md").write_text(
        "\n".join(report_lines), encoding="utf-8"
    )

    validation_columns = [
        "binary_cell_count_transfection_percent",
        "expression_brightness_weighted_efficiency_percent",
        "old_gfp_object_count_percent",
        "old_gfp_positive_area_percent",
        "mean_cell_gfp_score_0_to_255",
    ]
    validation_columns = [
        column for column in validation_columns if column in summary.columns
    ]
    correlations = summary[validation_columns].corr()

    def correlation(left: str, right: str) -> float:
        if left not in correlations.index or right not in correlations.columns:
            return float("nan")
        return float(correlations.loc[left, right])

    negative_summary = summary[
        summary["sample"].astype(str) == str(args.negative_control)
    ].iloc[0]
    validation = pd.DataFrame(
        [
            {
                "cells_classified": len(per_cell),
                "negative_control": str(args.negative_control),
                "negative_control_cells": len(negative),
                "negative_control_false_positive_percent": negative_summary[
                    "binary_cell_count_transfection_percent"
                ],
                "classification_threshold_0_to_255": threshold,
                "intensity_max_for_normalization": args.intensity_max,
                "new_count_correlation_with_mean_cell_gfp": correlation(
                    "binary_cell_count_transfection_percent",
                    "mean_cell_gfp_score_0_to_255",
                ),
                "old_object_count_correlation_with_mean_cell_gfp": correlation(
                    "old_gfp_object_count_percent",
                    "mean_cell_gfp_score_0_to_255",
                ),
                "old_area_correlation_with_mean_cell_gfp": correlation(
                    "old_gfp_positive_area_percent",
                    "mean_cell_gfp_score_0_to_255",
                ),
                "expression_efficiency_correlation_with_mean_cell_gfp": (
                    correlation(
                        "expression_brightness_weighted_efficiency_percent",
                        "mean_cell_gfp_score_0_to_255",
                    )
                ),
                "limitation": (
                    "Binary counting removes intensity weighting after classification, "
                    "but GFP must remain detectable above the negative-control distribution"
                ),
            }
        ]
    ).round(4)
    validation.to_csv(args.output_dir / "method_validation.csv", index=False)

    sensitivity_records: list[dict[str, float | int | str]] = []
    for quantile in (0.99, 0.995, 0.999):
        candidate_threshold = float(
            np.quantile(
                negative[score_column].to_numpy(), quantile, method="higher"
            )
        )
        for sample, group in per_cell.groupby("sample", sort=False):
            positives = int((group[score_column] > candidate_threshold).sum())
            expression_efficiency = float(
                np.clip(
                    (group[score_column] - candidate_threshold)
                    / (args.intensity_max - candidate_threshold),
                    0.0,
                    1.0,
                ).mean()
                * 100.0
            )
            sensitivity_records.append(
                {
                    "negative_control_quantile": quantile,
                    "classification_threshold_0_to_255": candidate_threshold,
                    "sample": sample,
                    "total_brightfield_cells": len(group),
                    "binary_gfp_positive_cells": positives,
                    "binary_cell_count_transfection_percent": positives
                    / len(group)
                    * 100.0,
                    "expression_brightness_weighted_efficiency_percent": (
                        expression_efficiency
                    ),
                }
            )
    pd.DataFrame(sensitivity_records).to_csv(
        args.output_dir / "threshold_sensitivity.csv", index=False
    )

    qc_samples = {str(args.negative_control)}
    ratios = summary.set_index("sample")["binary_cell_count_transfection_percent"]
    qc_samples.add(str((ratios - ratios.median()).abs().idxmin()))
    qc_samples.add(str(ratios.idxmax()))
    metadata = pd.DataFrame(image_metadata)
    for stale_qc in args.output_dir.glob("qc_*.png"):
        stale_qc.unlink()
    for sample in sorted(qc_samples):
        row = metadata[metadata["sample"].astype(str) == sample].iloc[0]
        sample_cells = per_cell[per_cell["sample"].astype(str) == sample]
        render_qc(
            args.input_dir / str(row["brightfield_tif"]),
            args.input_dir / str(row["gfp_tif"]),
            sample_cells,
            args.output_dir / f"qc_{sample}.png",
            args.crop_bottom,
        )

    print(f"Cells classified: {len(per_cell)}")
    print(f"Negative control: {args.negative_control} ({len(negative)} cells)")
    print(f"Threshold: {threshold:.4f} at quantile {args.negative_quantile}")
    print(f"Output: {args.output_dir.resolve()}")


if __name__ == "__main__":
    main()
