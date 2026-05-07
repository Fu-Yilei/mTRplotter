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
        "--use-fasta",
        action="store_true",
        help=(
            "Extract reads from trimmed_reads.fasta instead of trimmed_reads_to_poa.bam. "
            "This includes reads that failed to align to the POA consensus "
            "(e.g. rare somatic expansions that are much longer than the dominant allele). "
            "Applies to medaka samples only; TRGT samples are unaffected."
        ),
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
            "Use a moderate value on shared storage."
        ),
    )
    parser.add_argument(
        "--label-rotation",
        type=int,
        default=0,
        help="Rotation angle for x-axis tick labels in degrees. Default: 0 (horizontal).",
    )
    parser.add_argument(
        "--fontsize",
        type=int,
        default=11,
        help="Font size for x-axis tick labels. Axis labels and title scale proportionally. Default: 11.",
    )
    parser.add_argument(
        "--figure-width-per-sample",
        type=float,
        default=1.1,
        help="Figure width contribution per sample in inches. Default: 1.1.",
    )
    parser.add_argument(
        "--figure-height",
        type=float,
        default=8.0,
        help="Figure height in inches. Default: 8.0.",
    )
    parser.add_argument(
        "--no-haplotype-color",
        action="store_true",
        help="Plot all reads in a single color without haplotype distinction or legend.",
    )
    parser.add_argument(
        "--color-column",
        default=None,
        help=(
            "Color reads by a metadata column (e.g. donor). "
            "Overrides haplotype coloring. Each unique value gets a distinct color."
        ),
    )
    parser.add_argument(
        "--sort-columns",
        nargs="+",
        default=None,
        help="Sort samples on the x-axis by these metadata columns (e.g. donor source cell_type).",
    )
    parser.add_argument(
        "--sample-separators",
        action="store_true",
        help="Draw vertical separator lines between adjacent samples.",
    )
    parser.add_argument(
        "--separator-column",
        default=None,
        help=(
            "Draw vertical separator lines only when this sample metadata value "
            "changes between adjacent x-axis samples, e.g. donor."
        ),
    )
    parser.add_argument(
        "--group-label-column",
        default=None,
        help=(
            "Draw one centered label below each contiguous block where this "
            "sample metadata value is constant, e.g. donor."
        ),
    )
    parser.add_argument(
        "--subtract-reference-allele",
        action="store_true",
        help=(
            "Subtract the reference allele size (catalog end - start, looked up "
            "from --catalog-bed) from the y-axis read length on each plot. "
            "Requires --catalog-bed."
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
        use_fasta=args.use_fasta,
        label_rotation=args.label_rotation,
        fontsize=args.fontsize,
        no_haplotype_color=args.no_haplotype_color,
        color_column=args.color_column,
        sort_columns=args.sort_columns,
        sample_separators=args.sample_separators,
        separator_column=args.separator_column,
        group_label_column=args.group_label_column,
        figure_width_per_sample=args.figure_width_per_sample,
        figure_height=args.figure_height,
        subtract_reference_allele=args.subtract_reference_allele,
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
