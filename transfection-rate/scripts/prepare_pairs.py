#!/usr/bin/env python3
"""Create non-destructive, consistently named GFP/brightfield image pairs."""

from __future__ import annotations

import argparse
import csv
import os
import re
import shutil
from pathlib import Path


DEFAULT_GFP_REGEX = r"(?i)^.*?[-_](?P<sample>[^ _-]+)\s+(?:g|gfp|green)(?:[-_ ].*)?\.tiff?$"
DEFAULT_BRIGHTFIELD_REGEX = r"(?i)^.*?[-_](?P<sample>[^ _-]+)\s+(?:w|bf|brightfield|phase)(?:[-_ ].*)?\.tiff?$"


def parse_alias(value: str) -> tuple[str, str]:
    if "=" not in value:
        raise argparse.ArgumentTypeError("Aliases must use OLD=NEW syntax")
    old, new = value.split("=", 1)
    if not old.strip() or not new.strip():
        raise argparse.ArgumentTypeError("Alias labels cannot be empty")
    return old.strip().casefold(), new.strip()


def safe_sample_label(value: str) -> str:
    label = re.sub(r"[^A-Za-z0-9._-]+", "_", value.strip()).strip("._-")
    if not label:
        raise ValueError(f"Sample label {value!r} becomes empty after normalization")
    return label


def sample_sort_key(value: str) -> tuple[int, object, str]:
    if value.isdigit():
        return (0, int(value), value)
    return (1, value.casefold(), value)


def collect_matches(
    files: list[Path], pattern: re.Pattern[str], aliases: dict[str, str]
) -> dict[str, Path]:
    matches: dict[str, Path] = {}
    for path in files:
        match = pattern.match(path.name)
        if not match:
            continue
        raw_sample = match.groupdict().get("sample")
        if raw_sample is None:
            raise ValueError("Each channel regex must contain a named 'sample' group")
        sample = aliases.get(raw_sample.casefold(), raw_sample)
        if sample in matches:
            raise ValueError(
                f"Duplicate channel match for sample {sample!r}: "
                f"{matches[sample].name!r} and {path.name!r}"
            )
        matches[sample] = path
    return matches


def materialize(source: Path, destination: Path, mode: str) -> str:
    if destination.exists():
        if os.path.samefile(source, destination):
            return "existing-hardlink"
        raise FileExistsError(f"Refusing to overwrite {destination}")
    if mode == "copy":
        shutil.copy2(source, destination)
        return "copy"
    try:
        os.link(source, destination)
        return "hardlink"
    except OSError:
        shutil.copy2(source, destination)
        return "copy-fallback"


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Pair fluorescence and brightfield TIFF files without changing originals."
    )
    parser.add_argument("--input-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--gfp-regex", default=DEFAULT_GFP_REGEX)
    parser.add_argument("--brightfield-regex", default=DEFAULT_BRIGHTFIELD_REGEX)
    parser.add_argument(
        "--alias",
        action="append",
        default=[],
        type=parse_alias,
        metavar="OLD=NEW",
        help="Map an inconsistent source label to the intended sample label; repeat as needed.",
    )
    parser.add_argument(
        "--mode", choices=("hardlink", "copy"), default="hardlink"
    )
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()

    input_dir = args.input_dir.resolve()
    output_dir = args.output_dir.resolve()
    if not input_dir.is_dir():
        parser.error(f"Input directory does not exist: {input_dir}")
    if input_dir == output_dir:
        parser.error("Input and output directories must differ")

    aliases = dict(args.alias)
    try:
        gfp_pattern = re.compile(args.gfp_regex)
        brightfield_pattern = re.compile(args.brightfield_regex)
    except re.error as exc:
        parser.error(f"Invalid regular expression: {exc}")

    files = sorted(
        (path for path in input_dir.iterdir() if path.suffix.casefold() in {".tif", ".tiff"}),
        key=lambda path: path.name.casefold(),
    )
    gfp = collect_matches(files, gfp_pattern, aliases)
    brightfield = collect_matches(files, brightfield_pattern, aliases)

    missing_gfp = sorted(set(brightfield) - set(gfp), key=sample_sort_key)
    missing_brightfield = sorted(set(gfp) - set(brightfield), key=sample_sort_key)
    if missing_gfp or missing_brightfield:
        details = []
        if missing_gfp:
            details.append("missing fluorescence: " + ", ".join(missing_gfp))
        if missing_brightfield:
            details.append("missing brightfield: " + ", ".join(missing_brightfield))
        raise SystemExit("Pairing failed; " + "; ".join(details))
    if not gfp:
        raise SystemExit("No image pairs matched the supplied regular expressions")

    samples = sorted(gfp, key=sample_sort_key)
    width = max(2, len(str(len(samples) - 1)))
    rows: list[dict[str, str]] = []
    if not args.dry_run:
        output_dir.mkdir(parents=True, exist_ok=True)

    for index, sample in enumerate(samples):
        label = safe_sample_label(sample)
        gfp_name = f"{index:0{width}d}_{label}_g.tif"
        brightfield_name = f"{index:0{width}d}_{label}_w.tif"
        if args.dry_run:
            gfp_method = brightfield_method = "dry-run"
        else:
            gfp_method = materialize(gfp[sample], output_dir / gfp_name, args.mode)
            brightfield_method = materialize(
                brightfield[sample], output_dir / brightfield_name, args.mode
            )
        rows.append(
            {
                "sample": sample,
                "gfp_source": str(gfp[sample]),
                "brightfield_source": str(brightfield[sample]),
                "gfp_standardized": gfp_name,
                "brightfield_standardized": brightfield_name,
                "gfp_method": gfp_method,
                "brightfield_method": brightfield_method,
            }
        )

    manifest = output_dir / "pair_manifest.csv"
    if not args.dry_run:
        with manifest.open("w", newline="", encoding="utf-8-sig") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)

    print(f"Matched {len(rows)} complete image pairs")
    print(f"Manifest: {manifest if not args.dry_run else '(dry run; not written)'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
