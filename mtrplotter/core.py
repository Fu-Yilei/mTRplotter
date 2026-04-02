from __future__ import annotations

import csv
import gzip
import random
import re
import shutil
import statistics
import subprocess
from collections import defaultdict
from collections.abc import Iterator
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path
import tempfile

REGION_RE = re.compile(r"(chr[^_]+_\d+_\d+)")
HAP_RE = re.compile(r"_(hap\d+)_phased-set")
REGION_TEXT_RE = re.compile(r"^(chr[^:]+)[:_](\d+)[-_](\d+)$")
PLOT_PALETTE = {
    "hap1": "#4C78A8",
    "hap2": "#F58518",
    "unknown": "#888888",
}
REQUIRED_SAMPLE_COLUMNS = {"sample", "medaka_folder", "software", "flank_bp"}
CONTROL_SAMPLE_COLUMNS = {"sample", "medaka_folder", "software", "flank_bp"}
SUPPORTED_SOFTWARE = {"medaka", "trgt"}


def get_output_tmp_dir(output_dir: Path) -> Path:
    tmp_dir = output_dir / ".tmp"
    tmp_dir.mkdir(parents=True, exist_ok=True)
    return tmp_dir


@dataclass(frozen=True)
class TargetRegion:
    chrom: str
    start: int
    end: int
    region: str
    igv_locus: str

    @classmethod
    def from_region_text(cls, region_text: str) -> "TargetRegion":
        text = region_text.strip()
        match = REGION_TEXT_RE.fullmatch(text)
        if not match:
            raise ValueError(
                f"Could not parse region {region_text!r}. "
                "Expected chr:start-end or chr_start_end."
            )

        chrom, start_text, end_text = match.groups()
        start = int(start_text)
        end = int(end_text)
        if end <= start:
            raise ValueError(f"Region end must be larger than start: {region_text!r}")

        return cls(
            chrom=chrom,
            start=start,
            end=end,
            region=f"{chrom}_{start}_{end}",
            igv_locus=f"{chrom}:{start}-{end}",
        )

    @classmethod
    def from_bed_fields(cls, fields: list[str]) -> "TargetRegion":
        if len(fields) < 3:
            raise ValueError("BED rows must have at least three columns.")
        chrom = fields[0]
        start = int(fields[1])
        end = int(fields[2])
        return cls(
            chrom=chrom,
            start=start,
            end=end,
            region=f"{chrom}_{start}_{end}",
            igv_locus=f"{chrom}:{start}-{end}",
        )

    @property
    def slug(self) -> str:
        return self.region


@dataclass(frozen=True)
class SampleConfig:
    sample: str
    medaka_folder: Path
    software: str
    flank_bp: int
    metadata: dict[str, str]

    @property
    def ref_chunks_fasta(self) -> Path:
        return self.medaka_folder / "ref_chunks.fasta"

    @property
    def trimmed_reads_bam(self) -> Path:
        return self.medaka_folder / "trimmed_reads_to_poa.bam"

    @property
    def trimmed_reads_bam_index(self) -> Path:
        return self.medaka_folder / "trimmed_reads_to_poa.bam.bai"


@dataclass(frozen=True)
class TRGTSampleFiles:
    spanning_bam: Path
    vcf_gz: Path


def normalize_column_name(column_name: str) -> str:
    return column_name.strip().lower().replace(" ", "_")


def parse_non_negative_int(value: str, field_name: str, row_number: int) -> int:
    try:
        parsed_value = int(value)
    except ValueError as exc:
        raise ValueError(
            f"Sample table row {row_number} has a non-integer {field_name}: {value!r}"
        ) from exc

    if parsed_value < 0:
        raise ValueError(
            f"Sample table row {row_number} has a negative {field_name}: {parsed_value}"
        )
    return parsed_value


def deduplicate_paths(paths: list[Path]) -> list[Path]:
    deduplicated_paths: list[Path] = []
    seen_paths: set[str] = set()
    for path in paths:
        path_text = str(path)
        if path_text in seen_paths:
            continue
        seen_paths.add(path_text)
        deduplicated_paths.append(path)
    return deduplicated_paths


def select_existing_path(paths: list[Path]) -> Path | None:
    for path in paths:
        if path.is_file():
            return path
    return None


def resolve_trgt_sample_files(sample: SampleConfig) -> TRGTSampleFiles:
    candidate_bams: list[Path] = []
    candidate_vcfs: list[Path] = []
    sample_path = sample.medaka_folder

    if sample_path.is_dir():
        candidate_bams.extend(
            [
                sample_path / f"{sample.sample}.trgt.spanning.sorted.bam",
                sample_path / f"{sample.sample}.trgt.spanning.bam",
                *sorted(sample_path.glob(f"{sample.sample}*.trgt.spanning.sorted.bam")),
                *sorted(sample_path.glob(f"{sample.sample}*.trgt.spanning.bam")),
            ]
        )
        candidate_vcfs.extend(
            [
                sample_path / f"{sample.sample}.trgt.vcf.gz",
                *sorted(sample_path.glob(f"{sample.sample}*.trgt.vcf.gz")),
            ]
        )
    else:
        sample_path_text = str(sample_path)
        prefix_texts: list[str] = []
        if sample_path_text.endswith(".trgt.spanning.bam"):
            candidate_bams.append(Path(sample_path_text[: -len(".bam")] + ".sorted.bam"))
            candidate_bams.append(sample_path)
            prefix_texts.append(sample_path_text[: -len(".spanning.bam")])
        elif sample_path_text.endswith(".trgt.vcf.gz"):
            candidate_vcfs.append(sample_path)
            prefix_texts.append(sample_path_text[: -len(".vcf.gz")])
        elif sample_path_text.endswith(".trgt"):
            prefix_texts.append(sample_path_text)
        else:
            prefix_texts.extend([sample_path_text, f"{sample_path_text}.trgt"])

        for prefix_text in prefix_texts:
            candidate_bams.append(Path(f"{prefix_text}.spanning.sorted.bam"))
            candidate_bams.append(Path(f"{prefix_text}.spanning.bam"))
            candidate_vcfs.append(Path(f"{prefix_text}.vcf.gz"))

    bam_path = select_existing_path(deduplicate_paths(candidate_bams))
    if bam_path is None:
        raise FileNotFoundError(
            f"Could not resolve a TRGT spanning BAM for sample {sample.sample} from "
            f"{sample.medaka_folder}."
        )

    vcf_path = select_existing_path(deduplicate_paths(candidate_vcfs))
    if vcf_path is None:
        raise FileNotFoundError(
            f"Could not resolve a TRGT VCF.gz for sample {sample.sample} from "
            f"{sample.medaka_folder}."
        )

    return TRGTSampleFiles(spanning_bam=bam_path, vcf_gz=vcf_path)


def load_sample_table(
    sample_table_path: Path,
) -> list[SampleConfig]:
    samples: list[SampleConfig] = []

    with sample_table_path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"Sample table is empty: {sample_table_path}")

        normalized_names = {
            field_name: normalize_column_name(field_name)
            for field_name in reader.fieldnames
        }
        normalized_set = set(normalized_names.values())
        missing = REQUIRED_SAMPLE_COLUMNS - normalized_set
        if missing:
            raise ValueError(
                f"Sample table must contain columns {sorted(REQUIRED_SAMPLE_COLUMNS)}. "
                f"Missing: {sorted(missing)}"
            )

        for row_number, row in enumerate(reader, start=2):
            normalized_row = {
                normalized_names[key]: (value or "").strip()
                for key, value in row.items()
                if key is not None
            }
            sample_name = normalized_row["sample"]
            medaka_folder_text = normalized_row["medaka_folder"]
            if not sample_name or not medaka_folder_text:
                raise ValueError(
                    f"Sample table row {row_number} must include sample and medaka_folder."
                )

            software = (normalized_row.get("software") or "medaka").strip().lower()
            if software not in SUPPORTED_SOFTWARE:
                raise ValueError(
                    f"Sample table row {row_number} has unsupported software "
                    f"{software!r}. Expected one of {sorted(SUPPORTED_SOFTWARE)}."
                )

            flank_bp_text = normalized_row.get("flank_bp", "")
            if not flank_bp_text:
                raise ValueError(
                    f"Sample table row {row_number} must include flank_bp."
                )
            flank_bp = parse_non_negative_int(flank_bp_text, "flank_bp", row_number)

            merged_metadata = {
                key: value
                for key, value in normalized_row.items()
                if key not in CONTROL_SAMPLE_COLUMNS
            }

            medaka_folder = Path(medaka_folder_text)
            sample = SampleConfig(
                sample=sample_name,
                medaka_folder=medaka_folder,
                software=software,
                flank_bp=flank_bp,
                metadata=merged_metadata,
            )
            ensure_sample_paths_exist(sample)
            samples.append(sample)

    if not samples:
        raise ValueError(f"Sample table does not contain any samples: {sample_table_path}")
    return samples


def ensure_sample_paths_exist(sample: SampleConfig) -> None:
    if sample.software == "medaka":
        if not sample.medaka_folder.is_dir():
            raise FileNotFoundError(
                f"Medaka folder does not exist for sample {sample.sample}: {sample.medaka_folder}"
            )
        if not sample.ref_chunks_fasta.is_file():
            raise FileNotFoundError(
                f"Missing ref_chunks.fasta for sample {sample.sample}: {sample.ref_chunks_fasta}"
            )
        if not sample.trimmed_reads_bam.is_file():
            raise FileNotFoundError(
                f"Missing trimmed_reads_to_poa.bam for sample {sample.sample}: "
                f"{sample.trimmed_reads_bam}"
            )
        if not sample.trimmed_reads_bam_index.is_file():
            raise FileNotFoundError(
                f"Missing trimmed_reads_to_poa.bam.bai for sample {sample.sample}: "
                f"{sample.trimmed_reads_bam_index}"
            )
        return

    if sample.software == "trgt":
        resolve_trgt_sample_files(sample)
        return

    raise ValueError(
        f"Unsupported software {sample.software!r} for sample {sample.sample}."
    )


def load_targets(region_text: str | None, bed_path: Path | None) -> list[TargetRegion]:
    if (region_text is None) == (bed_path is None):
        raise ValueError("Provide exactly one of region_text or bed_path.")

    targets: list[TargetRegion] = []
    if region_text is not None:
        targets.append(TargetRegion.from_region_text(region_text))
    else:
        assert bed_path is not None
        with bed_path.open() as handle:
            for line_number, line in enumerate(handle, start=1):
                stripped = line.strip()
                if not stripped or stripped.startswith("#"):
                    continue
                fields = stripped.split("\t")
                try:
                    targets.append(TargetRegion.from_bed_fields(fields))
                except Exception as exc:
                    raise ValueError(
                        f"Could not parse BED line {line_number} in {bed_path}: {stripped}"
                    ) from exc

    unique_targets: dict[str, TargetRegion] = {}
    for target in targets:
        unique_targets[target.region] = target
    if not unique_targets:
        raise ValueError("No target loci were loaded.")
    return list(unique_targets.values())


def load_catalog_regions(catalog_bed_path: Path) -> set[str]:
    catalog_regions: set[str] = set()
    with catalog_bed_path.open() as handle:
        for line in handle:
            stripped = line.strip()
            if not stripped or stripped.startswith("#"):
                continue
            fields = stripped.split("\t")
            target = TargetRegion.from_bed_fields(fields)
            catalog_regions.add(target.region)
    return catalog_regions


def resolve_jobs(jobs: int, item_count: int) -> int:
    if jobs < 1:
        raise ValueError(f"jobs must be at least 1, got {jobs}")
    return min(jobs, max(1, item_count))


def adjust_read_length(raw_read_length: int, flank_bp: int) -> int:
    if flank_bp < 0:
        raise ValueError(f"flank_bp must be non-negative, got {flank_bp}")
    return raw_read_length - (2 * flank_bp)


def ensure_bam_tools_available() -> None:
    if shutil.which("samtools") is None:
        raise RuntimeError("samtools is required because mTRplotter extracts reads from BAM.")


def find_target_chunks_in_ref_chunks(
    ref_chunks_path: Path,
    target_regions: set[str],
) -> dict[str, list[str]]:
    found_chunks: dict[str, list[str]] = defaultdict(list)
    with ref_chunks_path.open() as handle:
        for line in handle:
            if not line.startswith(">"):
                continue
            match = REGION_RE.search(line)
            if not match:
                continue
            region = match.group(1)
            if region in target_regions:
                found_chunks[region].append(line[1:].strip())
    return dict(found_chunks)


def scan_sample_ref_chunks(
    sample: SampleConfig,
    target_regions: set[str],
) -> tuple[str, dict[str, list[str]]]:
    found_chunks = find_target_chunks_in_ref_chunks(sample.ref_chunks_fasta, target_regions)
    return sample.sample, found_chunks


def parse_vcf_info_field(info_text: str) -> dict[str, str]:
    info_map: dict[str, str] = {}
    for entry in info_text.split(";"):
        if not entry:
            continue
        if "=" in entry:
            key, value = entry.split("=", 1)
            info_map[key] = value
        else:
            info_map[entry] = ""
    return info_map


def find_target_trids_in_vcf(
    vcf_gz_path: Path,
    target_regions: set[str],
) -> dict[str, list[str]]:
    found_trids: dict[str, list[str]] = defaultdict(list)
    with gzip.open(vcf_gz_path, "rt") as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 8:
                continue

            chrom = fields[0]
            start = int(fields[1])
            info_map = parse_vcf_info_field(fields[7])
            end_text = info_map.get("END")
            trid = info_map.get("TRID")
            if not end_text or not trid:
                continue

            region = f"{chrom}_{start}_{int(end_text)}"
            if region not in target_regions:
                continue
            if trid not in found_trids[region]:
                found_trids[region].append(trid)

    return dict(found_trids)


def scan_sample_targets(
    sample: SampleConfig,
    target_regions: set[str],
) -> tuple[str, dict[str, list[str]]]:
    if sample.software == "medaka":
        return scan_sample_ref_chunks(sample, target_regions)

    trgt_files = resolve_trgt_sample_files(sample)
    return sample.sample, find_target_trids_in_vcf(trgt_files.vcf_gz, target_regions)


def validate_targets(
    samples: list[SampleConfig],
    targets: list[TargetRegion],
    catalog_bed_path: Path | None,
    require_region_in_all_samples: bool,
    jobs: int,
) -> tuple[list[dict[str, str]], dict[str, dict[str, list[str]]]]:
    target_regions = {target.region for target in targets}
    validation_rows: list[dict[str, str]] = []

    if catalog_bed_path is not None:
        catalog_regions = load_catalog_regions(catalog_bed_path)
        missing_from_catalog = sorted(target_regions - catalog_regions)
        if missing_from_catalog:
            raise ValueError(
                "Requested loci are not a subset of the provided catalog BED: "
                + ", ".join(missing_from_catalog)
            )

    found_chunks_by_sample: dict[str, dict[str, list[str]]] = {}
    worker_count = resolve_jobs(jobs, len(samples))

    if worker_count == 1:
        for sample in samples:
            sample_name, found_entries = scan_sample_targets(sample, target_regions)
            found_chunks_by_sample[sample_name] = found_entries
    else:
        with ThreadPoolExecutor(max_workers=worker_count) as executor:
            futures = {
                executor.submit(scan_sample_targets, sample, target_regions): sample.sample
                for sample in samples
            }
            for completed_index, future in enumerate(as_completed(futures), start=1):
                sample_name, found_entries = future.result()
                found_chunks_by_sample[sample_name] = found_entries
                print(
                    f"[validate] completed {completed_index}/{len(samples)}: {sample_name}",
                    flush=True,
                )

    missing_by_sample: dict[str, list[str]] = {}
    for sample in samples:
        found_chunks = found_chunks_by_sample[sample.sample]
        found_regions = set(found_chunks)
        missing_regions = sorted(target_regions - found_regions)
        missing_by_sample[sample.sample] = missing_regions
        for target in targets:
            validation_rows.append(
                {
                    "sample": sample.sample,
                    "software": sample.software,
                    "source_path": str(sample.medaka_folder),
                    "region": target.region,
                    "in_source": str(target.region in found_regions).lower(),
                    "matching_entries": ",".join(found_chunks.get(target.region, [])),
                }
            )

    if require_region_in_all_samples:
        missing_messages = []
        for sample_name, missing_regions in missing_by_sample.items():
            if missing_regions:
                missing_messages.append(f"{sample_name}: {', '.join(missing_regions)}")
        if missing_messages:
            raise ValueError(
                "Requested loci are absent from some sample inputs: "
                + "; ".join(missing_messages)
            )

    return validation_rows, found_chunks_by_sample


def parse_haplotype(text: str) -> str:
    match = HAP_RE.search(text)
    if match:
        return match.group(1)
    return "unknown"


def parse_optional_tags(optional_fields: list[str]) -> dict[str, str]:
    tags: dict[str, str] = {}
    for field in optional_fields:
        parts = field.split(":", 2)
        if len(parts) != 3:
            continue
        tags[parts[0]] = parts[2]
    return tags


def stream_command_lines(command: list[str]) -> Iterator[str]:
    process = subprocess.Popen(
        command,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )
    assert process.stdout is not None
    assert process.stderr is not None

    try:
        for line in process.stdout:
            yield line.rstrip("\n")
    finally:
        process.stdout.close()
        try:
            stderr_text = process.stderr.read()
        finally:
            process.stderr.close()
        return_code = process.wait()
        if return_code != 0:
            raise subprocess.CalledProcessError(
                returncode=return_code,
                cmd=command,
                stderr=stderr_text,
            )


def parse_trgt_haplotype(tags: dict[str, str]) -> str:
    hap_value = tags.get("HP", "")
    if hap_value == "1":
        return "hap1"
    if hap_value == "2":
        return "hap2"
    return "unknown"


def parse_trgt_flanks(
    fl_tag_value: str | None,
    fallback_flank_bp: int,
) -> tuple[int, int]:
    if fl_tag_value:
        parts = fl_tag_value.split(",")
        value_parts = parts[1:] if parts and parts[0].isalpha() else parts
        if len(value_parts) >= 2:
            try:
                return int(value_parts[0]), int(value_parts[1])
            except ValueError:
                pass
    return fallback_flank_bp, fallback_flank_bp


def collect_sample_read_rows_from_bam(
    sample: SampleConfig,
    target_lookup: dict[str, TargetRegion],
    target_chunks_by_region: dict[str, list[str]],
    keep_length_one: bool,
) -> tuple[str, list[dict[str, str | int]]]:
    chunk_names = [
        chunk_name
        for region in target_lookup
        for chunk_name in target_chunks_by_region.get(region, [])
    ]
    if not chunk_names:
        return sample.sample, []

    command = [
        "samtools",
        "view",
        "-F",
        "2304",
        str(sample.trimmed_reads_bam),
        *chunk_names,
    ]
    read_rows: list[dict[str, str | int]] = []
    for line in stream_command_lines(command):
        if not line:
            continue
        fields = line.split("\t")
        if len(fields) < 11:
            continue

        reference_name = fields[2]
        sequence = fields[9]
        if sequence == "*":
            continue

        region_match = REGION_RE.search(reference_name)
        if not region_match:
            continue

        region = region_match.group(1)
        target = target_lookup[region]
        raw_read_length = len(sequence)
        if not keep_length_one and raw_read_length == 1:
            continue
        adjusted_read_length = adjust_read_length(raw_read_length, sample.flank_bp)

        row: dict[str, str | int] = {
            "sample": sample.sample,
            "medaka_folder": str(sample.medaka_folder),
            "software": sample.software,
            "chrom": target.chrom,
            "start": target.start,
            "end": target.end,
            "region": target.region,
            "igv_locus": target.igv_locus,
            "hap": parse_haplotype(reference_name),
            "raw_read_length": raw_read_length,
            "flank_bp": sample.flank_bp,
            "flank_left_bp": sample.flank_bp,
            "flank_right_bp": sample.flank_bp,
            "read_length": adjusted_read_length,
        }
        for key, value in sample.metadata.items():
            row[key] = value
        read_rows.append(row)

    return sample.sample, read_rows


def collect_sample_read_rows_from_trgt_bam(
    sample: SampleConfig,
    target_lookup: dict[str, TargetRegion],
    target_trids_by_region: dict[str, list[str]],
    keep_length_one: bool,
    temp_dir: Path,
) -> tuple[str, list[dict[str, str | int]]]:
    trgt_files = resolve_trgt_sample_files(sample)

    region_by_trid: dict[str, str] = {}
    for region, trids in target_trids_by_region.items():
        unique_trids = sorted(set(trids))
        if len(unique_trids) > 1:
            raise ValueError(
                f"TRGT sample {sample.sample} matched multiple TRIDs for region {region}: "
                f"{', '.join(unique_trids)}"
            )
        if unique_trids:
            region_by_trid[unique_trids[0]] = region

    if not region_by_trid:
        return sample.sample, []

    command: list[str]
    temporary_bed_path: Path | None = None
    bam_index_path = Path(f"{trgt_files.spanning_bam}.bai")
    if not bam_index_path.is_file():
        print(
            f"[warning] TRGT spanning BAM for {sample.sample} has no index "
            f"({trgt_files.spanning_bam}); sort and index it with "
            f"'samtools sort' + 'samtools index' to avoid scanning the full BAM file.",
            flush=True,
        )
    if bam_index_path.is_file():
        with tempfile.NamedTemporaryFile(
            mode="w",
            suffix=".bed",
            prefix=f"{sample.sample}.",
            dir=temp_dir,
            delete=False,
        ) as temporary_bed_handle:
            for region in sorted(set(region_by_trid.values())):
                target = target_lookup[region]
                temporary_bed_handle.write(
                    f"{target.chrom}\t{target.start}\t{target.end}\n"
                )
        temporary_bed_path = Path(temporary_bed_handle.name)
        command = [
            "samtools",
            "view",
            "-F",
            "2304",
            "-M",
            "-L",
            str(temporary_bed_path),
            str(trgt_files.spanning_bam),
        ]
    elif len(region_by_trid) == 1:
        trid_filter = next(iter(region_by_trid))
        command = [
            "samtools",
            "view",
            "-F",
            "2304",
            "-d",
            f"TR:{trid_filter}",
            str(trgt_files.spanning_bam),
        ]
    else:
        command = [
            "samtools",
            "view",
            "-F",
            "2304",
            str(trgt_files.spanning_bam),
        ]

    read_rows: list[dict[str, str | int]] = []
    try:
        for line in stream_command_lines(command):
            if not line:
                continue
            fields = line.split("\t")
            if len(fields) < 11:
                continue

            sequence = fields[9]
            if sequence == "*":
                continue

            tags = parse_optional_tags(fields[11:])
            trid = tags.get("TR")
            if not trid:
                continue

            region = region_by_trid.get(trid)
            if region is None:
                continue

            raw_read_length = len(sequence)
            if not keep_length_one and raw_read_length == 1:
                continue

            flank_left_bp, flank_right_bp = parse_trgt_flanks(
                tags.get("FL"),
                sample.flank_bp,
            )
            adjusted_read_length = raw_read_length - flank_left_bp - flank_right_bp
            target = target_lookup[region]
            representative_flank_bp = flank_left_bp if flank_left_bp == flank_right_bp else ""

            row: dict[str, str | int] = {
                "sample": sample.sample,
                "medaka_folder": str(sample.medaka_folder),
                "software": sample.software,
                "chrom": target.chrom,
                "start": target.start,
                "end": target.end,
                "region": target.region,
                "igv_locus": target.igv_locus,
                "hap": parse_trgt_haplotype(tags),
                "raw_read_length": raw_read_length,
                "flank_bp": representative_flank_bp,
                "flank_left_bp": flank_left_bp,
                "flank_right_bp": flank_right_bp,
                "read_length": adjusted_read_length,
            }
            for key, value in sample.metadata.items():
                row[key] = value
            read_rows.append(row)
    finally:
        if temporary_bed_path is not None and temporary_bed_path.exists():
            temporary_bed_path.unlink()

    return sample.sample, read_rows


def collect_rows_for_sample(
    sample: SampleConfig,
    target_lookup: dict[str, TargetRegion],
    keep_length_one: bool,
    target_matches_by_region: dict[str, list[str]],
    temp_dir: Path,
) -> tuple[str, list[dict[str, str | int]]]:
    if sample.software == "trgt":
        return collect_sample_read_rows_from_trgt_bam(
            sample=sample,
            target_lookup=target_lookup,
            target_trids_by_region=target_matches_by_region,
            keep_length_one=keep_length_one,
            temp_dir=temp_dir,
        )

    return collect_sample_read_rows_from_bam(
        sample=sample,
        target_lookup=target_lookup,
        target_chunks_by_region=target_matches_by_region,
        keep_length_one=keep_length_one,
    )


def collect_read_rows(
    samples: list[SampleConfig],
    targets: list[TargetRegion],
    keep_length_one: bool,
    jobs: int,
    target_chunks_by_sample: dict[str, dict[str, list[str]]],
    temp_dir: Path,
) -> list[dict[str, str | int]]:
    target_lookup = {target.region: target for target in targets}
    read_rows: list[dict[str, str | int]] = []
    worker_count = resolve_jobs(jobs, len(samples))

    sample_rows: dict[str, list[dict[str, str | int]]] = {}
    if worker_count == 1:
        for sample in samples:
            sample_name, rows = collect_rows_for_sample(
                sample=sample,
                target_lookup=target_lookup,
                keep_length_one=keep_length_one,
                target_matches_by_region=target_chunks_by_sample[sample.sample],
                temp_dir=temp_dir,
            )
            sample_rows[sample_name] = rows
    else:
        with ThreadPoolExecutor(max_workers=worker_count) as executor:
            futures = {}
            for sample in samples:
                future = executor.submit(
                    collect_rows_for_sample,
                    sample,
                    target_lookup,
                    keep_length_one,
                    target_chunks_by_sample[sample.sample],
                    temp_dir,
                )
                futures[future] = sample.sample
            for completed_index, future in enumerate(as_completed(futures), start=1):
                sample_name, rows = future.result()
                sample_rows[sample_name] = rows
                print(
                    f"[collect] completed {completed_index}/{len(samples)}: {sample_name}",
                    flush=True,
                )

    for sample in samples:
        read_rows.extend(sample_rows[sample.sample])

    return read_rows


def summarize_read_rows(read_rows: list[dict[str, str | int]]) -> list[dict[str, str | int | float]]:
    grouped_lengths: dict[tuple[str, str, str], list[int]] = defaultdict(list)
    row_lookup: dict[tuple[str, str, str], dict[str, str | int]] = {}

    for row in read_rows:
        key = (str(row["region"]), str(row["sample"]), str(row["hap"]))
        grouped_lengths[key].append(int(row["read_length"]))
        row_lookup[key] = row

    summary_rows: list[dict[str, str | int | float]] = []
    for key in sorted(grouped_lengths):
        lengths = grouped_lengths[key]
        example_row = row_lookup[key]
        summary_row: dict[str, str | int | float] = {
            "sample": example_row["sample"],
            "medaka_folder": example_row["medaka_folder"],
            "chrom": example_row["chrom"],
            "start": example_row["start"],
            "end": example_row["end"],
            "region": example_row["region"],
            "igv_locus": example_row["igv_locus"],
            "hap": example_row["hap"],
            "flank_bp": example_row["flank_bp"],
            "n_reads": len(lengths),
            "mean_length": round(statistics.fmean(lengths), 3),
            "median_length": statistics.median(lengths),
            "min_length": min(lengths),
            "max_length": max(lengths),
        }
        if len(lengths) > 1:
            summary_row["stdev_length"] = round(statistics.stdev(lengths), 3)
        else:
            summary_row["stdev_length"] = ""
        for field_name, field_value in example_row.items():
            if field_name in summary_row or field_name in {"read_length", "raw_read_length"}:
                continue
            summary_row[field_name] = field_value
        summary_rows.append(summary_row)

    return summary_rows


def build_sample_labels(samples: list[SampleConfig], label_columns: list[str]) -> dict[str, str]:
    labels: dict[str, str] = {}
    normalized_label_columns = [normalize_column_name(column) for column in label_columns]

    available_columns = {"sample", "software", "flank_bp"}
    for sample in samples:
        available_columns.update(sample.metadata.keys())
    missing_columns = sorted(column for column in normalized_label_columns if column not in available_columns)
    if missing_columns:
        raise ValueError(
            "Requested label columns are not available after metadata merge: "
            + ", ".join(missing_columns)
        )

    for sample in samples:
        label_parts: list[str] = []
        for column in normalized_label_columns:
            if column == "sample":
                label_parts.append(sample.sample)
            elif column == "software":
                label_parts.append(sample.software)
            elif column == "flank_bp":
                label_parts.append(str(sample.flank_bp))
            else:
                label_parts.append(sample.metadata.get(column, ""))

        label_text = "\n".join(part for part in label_parts if part)
        labels[sample.sample] = label_text or sample.sample
    return labels


def describe_y_axis(region_rows: list[dict[str, str | int]]) -> str:
    flank_pairs = {
        (int(row.get("flank_left_bp", 0)), int(row.get("flank_right_bp", 0)))
        for row in region_rows
    }
    if not flank_pairs or flank_pairs == {(0, 0)}:
        return "Read length"

    if len(flank_pairs) == 1:
        flank_left_bp, flank_right_bp = next(iter(flank_pairs))
        return (
            "Read length minus "
            f"{flank_left_bp + flank_right_bp} bp flanks"
        )

    return "Flank-adjusted read length"


def plot_target_region(
    read_rows: list[dict[str, str | int]],
    target: TargetRegion,
    samples: list[SampleConfig],
    output_dir: Path,
    label_columns: list[str],
    figure_prefix: str,
    random_seed: int = 42,
) -> Path:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    sample_order = [sample.sample for sample in samples]
    sample_positions = {sample_name: index for index, sample_name in enumerate(sample_order)}
    sample_labels = build_sample_labels(samples, label_columns)
    region_rows = [row for row in read_rows if row["region"] == target.region]

    figure_width = max(12.0, len(samples) * 0.7)
    fig, ax = plt.subplots(figsize=(figure_width, 7))
    rng = random.Random(f"{random_seed}:{target.region}")

    for hap in ("hap1", "hap2", "unknown"):
        xs: list[float] = []
        ys: list[int] = []
        for row in region_rows:
            if row["hap"] != hap:
                continue
            base_x = sample_positions[str(row["sample"])]
            xs.append(base_x + rng.uniform(-0.22, 0.22))
            ys.append(int(row["read_length"]))
        if xs:
            ax.scatter(
                xs,
                ys,
                s=34,
                alpha=0.75,
                c=PLOT_PALETTE[hap],
                edgecolors="black",
                linewidths=0.2,
                label=hap,
            )

    ax.set_xticks(list(range(len(sample_order))))
    ax.set_xticklabels(
        [sample_labels[sample_name] for sample_name in sample_order],
        rotation=90,
        ha="center",
        fontsize=8,
    )
    ax.set_xlabel("Sample")
    ax.set_ylabel(describe_y_axis(region_rows))
    ax.set_title(f"Read length distribution across samples\n{target.igv_locus}")
    ax.grid(axis="y", color="#D9D9D9", linewidth=0.7)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    handles, labels = ax.get_legend_handles_labels()
    if handles:
        unique_labels: dict[str, object] = {}
        for handle, label in zip(handles, labels):
            unique_labels[label] = handle
        ax.legend(
            unique_labels.values(),
            unique_labels.keys(),
            title="Haplotype",
            bbox_to_anchor=(1.01, 1),
            loc="upper left",
            frameon=True,
        )

    plt.tight_layout()
    output_path = output_dir / f"{figure_prefix}.{target.slug}.png"
    fig.savefig(output_path, bbox_inches="tight")
    plt.close(fig)
    return output_path


def write_tsv(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        path.write_text("")
        return

    fieldnames = list(rows[0].keys())
    opener = gzip.open if path.suffix == ".gz" else open
    open_mode = "wt" if path.suffix == ".gz" else "w"
    with opener(path, open_mode, newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def run_workflow(
    sample_table_path: Path,
    output_dir: Path,
    region_text: str | None = None,
    bed_path: Path | None = None,
    catalog_bed_path: Path | None = None,
    label_columns: list[str] | None = None,
    keep_length_one: bool = False,
    require_region_in_all_samples: bool = True,
    figure_prefix: str = "read_length_by_sample",
    jobs: int = 1,
) -> dict[str, object]:
    labels = label_columns or ["sample"]
    ensure_bam_tools_available()
    samples = load_sample_table(sample_table_path=sample_table_path)
    targets = load_targets(region_text=region_text, bed_path=bed_path)
    output_dir.mkdir(parents=True, exist_ok=True)
    temp_dir = get_output_tmp_dir(output_dir)

    validation_rows, target_chunks_by_sample = validate_targets(
        samples=samples,
        targets=targets,
        catalog_bed_path=catalog_bed_path,
        require_region_in_all_samples=require_region_in_all_samples,
        jobs=jobs,
    )
    read_rows = collect_read_rows(
        samples=samples,
        targets=targets,
        keep_length_one=keep_length_one,
        jobs=jobs,
        target_chunks_by_sample=target_chunks_by_sample,
        temp_dir=temp_dir,
    )
    if not read_rows:
        raise ValueError("No reads matched the requested loci.")

    summary_rows = summarize_read_rows(read_rows)
    figure_paths = [
        plot_target_region(
            read_rows=read_rows,
            target=target,
            samples=samples,
            output_dir=output_dir,
            label_columns=labels,
            figure_prefix=figure_prefix,
        )
        for target in targets
    ]

    reads_tsv_path = output_dir / "per_read_lengths.tsv.gz"
    summary_tsv_path = output_dir / "read_length_summary.tsv"
    validation_tsv_path = output_dir / "region_validation.tsv"
    write_tsv(reads_tsv_path, read_rows)
    write_tsv(summary_tsv_path, summary_rows)
    write_tsv(validation_tsv_path, validation_rows)

    return {
        "figure_paths": figure_paths,
        "reads_tsv_path": reads_tsv_path,
        "summary_tsv_path": summary_tsv_path,
        "validation_tsv_path": validation_tsv_path,
        "read_source": "bam",
    }
