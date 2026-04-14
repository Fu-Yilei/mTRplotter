from __future__ import annotations

import argparse
from pathlib import Path

from .core import run_workflow


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Plot locus-specific read-length distributions across samples from "
            "medaka or TRGT outputs."
        )
    )
    parser.add_argument(
        "--sample-table",
        required=True,
        type=Path,
        help=(
            "Tab-delimited table with required columns sample, medaka_folder, "
            "software, flank_bp, and any additional metadata or label columns."
        ),
    )

    target_group = parser.add_mutually_exclusive_group(required=True)
    target_group.add_argument(
        "--region",
        help="Single target locus, for example chr1:99682978-99683197.",
    )
    target_group.add_argument(
        "--bed",
        type=Path,
        help="BED file containing target loci. Only the first three columns are used.",
    )

    parser.add_argument(
        "--output-dir",
        type=Path,
        help="Directory for figures, tables, params.json, and log.out. Default: current working directory.",
    )
    parser.add_argument(
        "--catalog-bed",
        type=Path,
        help="Optional BED file for an explicit catalog subset check.",
    )
    parser.add_argument(
        "--label-columns",
        nargs="+",
        default=["sample"],
        help="Columns to show in the x-axis labels. Default: sample",
    )
    parser.add_argument(
        "--keep-length-one",
        action="store_true",
        help="Keep reads with length 1 instead of filtering them out.",
    )
    parser.add_argument(
        "--allow-missing-regions-in-samples",
        action="store_true",
        help="Allow requested loci to be absent from some sample inputs.",
    )
    parser.add_argument(
        "--figure-prefix",
        default="read_length_by_sample",
        help="Filename prefix for saved figures.",
    )
    parser.add_argument(
        "--jobs",
        type=int,
        default=1,
        help=(
            "Number of samples to process in parallel. "
            "Use a moderate value on shared storage because all extraction runs through BAMs."
        ),
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    output_dir = (args.output_dir or Path.cwd()).resolve()

    result = run_workflow(
        sample_table_path=args.sample_table,
        output_dir=output_dir,
        region_text=args.region,
        bed_path=args.bed,
        catalog_bed_path=args.catalog_bed,
        label_columns=args.label_columns,
        keep_length_one=args.keep_length_one,
        require_region_in_all_samples=not args.allow_missing_regions_in_samples,
        figure_prefix=args.figure_prefix,
        jobs=args.jobs,
    )

    print(f"Figures written: {len(result['figure_paths'])}")
    for figure_path in result["figure_paths"]:
        print(figure_path)
    print(f"Per-read table: {result['reads_tsv_path']}")
    print(f"Summary table: {result['summary_tsv_path']}")
    print(f"Validation table: {result['validation_tsv_path']}")
    print(f"Parameters: {result['params_path']}")
    print(f"Log: {result['log_path']}")
    return 0
