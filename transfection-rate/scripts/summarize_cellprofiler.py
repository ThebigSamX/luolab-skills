#!/usr/bin/env python3
"""Create a clean transfection-rate summary from CellProfiler Image.csv."""

from __future__ import annotations

import argparse
import csv
import math
import re
from pathlib import Path


REQUIRED_COLUMNS = {
    "Count_PositiveCells",
    "Count_TotalCells",
    "AreaOccupied_AreaOccupied_GFPPositive",
    "AreaOccupied_TotalArea_GFPPositive",
}
INTERPRETATION = (
    "Estimated from brightfield cell centers without nuclear stain; suitable for "
    "relative comparison, not a validated absolute transfection efficiency"
)


def read_csv(path: Path) -> tuple[list[dict[str, str]], str]:
    raw = path.read_bytes()
    last_error: UnicodeDecodeError | None = None
    for encoding in ("utf-8-sig", "utf-8", "gb18030"):
        try:
            text = raw.decode(encoding)
        except UnicodeDecodeError as exc:
            last_error = exc
            continue
        rows = list(csv.DictReader(text.splitlines()))
        if rows:
            return rows, encoding
    if last_error:
        raise last_error
    raise ValueError(f"No data rows found in {path}")


def number(row: dict[str, str], column: str) -> float:
    value = row.get(column, "").strip()
    try:
        result = float(value)
    except ValueError as exc:
        raise ValueError(f"Invalid numeric value in {column}: {value!r}") from exc
    if not math.isfinite(result):
        raise ValueError(f"Non-finite numeric value in {column}: {value!r}")
    return result


def integer(row: dict[str, str], column: str) -> int:
    value = number(row, column)
    rounded = round(value)
    if not math.isclose(value, rounded, abs_tol=1e-6):
        raise ValueError(f"Expected integer count in {column}, got {value}")
    return int(rounded)


def sample_name(row: dict[str, str]) -> str:
    metadata = row.get("Metadata_Sample", "").strip()
    if metadata:
        return metadata
    filename = row.get("FileName_GFPColor", "").strip()
    match = re.match(r"^\d+_(?P<sample>.+)_g\.tiff?$", filename, flags=re.IGNORECASE)
    if match:
        return match.group("sample")
    return Path(filename).stem or row.get("ImageNumber", "unknown")


def first_present_number(row: dict[str, str], columns: tuple[str, ...]) -> float | None:
    for column in columns:
        if row.get(column, "").strip():
            return number(row, column)
    return None


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Summarize CellProfiler fluorescence area and estimated cell-count ratios."
    )
    parser.add_argument("--image-csv", type=Path, required=True)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    image_csv = args.image_csv.resolve()
    if not image_csv.is_file():
        parser.error(f"Image.csv does not exist: {image_csv}")
    output = (
        args.output.resolve()
        if args.output
        else image_csv.with_name("transfection_area_and_cell_count_summary.csv")
    )

    rows, encoding = read_csv(image_csv)
    missing = sorted(REQUIRED_COLUMNS - set(rows[0]))
    if missing:
        raise SystemExit("Missing required CellProfiler columns: " + ", ".join(missing))

    summaries: list[dict[str, object]] = []
    violations: list[str] = []
    total_module_errors = 0
    for row in rows:
        sample = sample_name(row)
        total = integer(row, "Count_TotalCells")
        positive = integer(row, "Count_PositiveCells")
        occupied = number(row, "AreaOccupied_AreaOccupied_GFPPositive")
        analyzed_area = number(row, "AreaOccupied_TotalArea_GFPPositive")
        if total <= 0:
            raise SystemExit(f"Sample {sample}: total-cell count must be positive")
        if analyzed_area <= 0:
            raise SystemExit(f"Sample {sample}: analyzed image area must be positive")
        if positive > total:
            violations.append(sample)

        module_errors = sum(
            integer(row, column)
            for column in row
            if column.startswith("ModuleError_") and row[column].strip()
        )
        total_module_errors += module_errors
        threshold = first_present_number(
            row,
            (
                "Threshold_FinalThreshold_GFPPositive",
                "Threshold_OrigThreshold_GFPPositive",
                "Threshold_FinalThreshold_PositiveCells",
            ),
        )
        summaries.append(
            {
                "sample": sample,
                "estimated_total_cells": total,
                "estimated_gfp_positive_cells": positive,
                "estimated_cell_count_ratio_percent": positive / total * 100.0,
                "gfp_positive_area_percent": occupied / analyzed_area * 100.0,
                "fixed_gfp_threshold_0_to_1": threshold,
                "module_error_count": module_errors,
                "paired_gfp_tif": row.get("FileName_GFPColor", ""),
                "paired_brightfield_tif": row.get("FileName_Brightfield", ""),
            }
        )

    ranked = sorted(
        summaries,
        key=lambda item: (-float(item["estimated_cell_count_ratio_percent"]), str(item["sample"])),
    )
    ranks = {id(item): rank for rank, item in enumerate(ranked, start=1)}
    output_rows: list[dict[str, object]] = []
    for item in summaries:
        output_rows.append(
            {
                "sample": item["sample"],
                "rank_by_cell_count_ratio": ranks[id(item)],
                "estimated_total_cells": item["estimated_total_cells"],
                "estimated_gfp_positive_cells": item["estimated_gfp_positive_cells"],
                "estimated_cell_count_ratio_percent": f'{float(item["estimated_cell_count_ratio_percent"]):.3f}',
                "gfp_positive_area_percent": f'{float(item["gfp_positive_area_percent"]):.3f}',
                "fixed_gfp_threshold_0_to_1": (
                    "" if item["fixed_gfp_threshold_0_to_1"] is None else f'{float(item["fixed_gfp_threshold_0_to_1"]):.4f}'
                ),
                "module_error_count": item["module_error_count"],
                "paired_gfp_tif": item["paired_gfp_tif"],
                "paired_brightfield_tif": item["paired_brightfield_tif"],
                "interpretation": INTERPRETATION,
            }
        )

    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(output_rows[0]))
        writer.writeheader()
        writer.writerows(output_rows)

    print(f"Read {len(rows)} image sets using {encoding}")
    print(f"Module errors: {total_module_errors}")
    print(f"Positive count exceeds total count: {len(violations)}")
    if violations:
        print("Flagged samples: " + ", ".join(violations))
    print(f"Summary: {output}")
    return 2 if total_module_errors or violations else 0


if __name__ == "__main__":
    raise SystemExit(main())
