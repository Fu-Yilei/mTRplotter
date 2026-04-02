from __future__ import annotations

import csv
import gzip
import shutil
import subprocess
import tempfile
import unittest
from pathlib import Path

from mtrplotter.core import TargetRegion, run_workflow


class WorkflowTests(unittest.TestCase):
    def setUp(self) -> None:
        self.tempdir = tempfile.TemporaryDirectory()
        self.root = Path(self.tempdir.name)

    def tearDown(self) -> None:
        self.tempdir.cleanup()

    def write_text(self, path: Path, text: str) -> None:
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(text)

    def write_gzip_text(self, path: Path, text: str) -> None:
        path.parent.mkdir(parents=True, exist_ok=True)
        with gzip.open(path, "wt") as handle:
            handle.write(text)

    def make_medaka_bam_sample(
        self,
        sample_name: str,
        ref_headers: list[str],
        sam_records: list[str],
    ) -> Path:
        medaka_dir = self.root / f"{sample_name}.medaka"
        self.write_text(medaka_dir / "ref_chunks.fasta", "".join(ref_headers))
        sam_header = "@HD\tVN:1.6\tSO:coordinate\n" + "".join(
            f"@SQ\tSN:{header[1:].strip()}\tLN:1000\n" for header in ref_headers
        )
        sam_text = sam_header + "".join(sam_records)
        sam_path = medaka_dir / "trimmed_reads_to_poa.sam"
        bam_path = medaka_dir / "trimmed_reads_to_poa.bam"
        self.write_text(sam_path, sam_text)
        subprocess.run(
            [
                "samtools",
                "view",
                "-b",
                "-o",
                str(bam_path),
                str(sam_path),
            ],
            check=True,
        )
        subprocess.run(["samtools", "index", str(bam_path)], check=True)
        return medaka_dir

    def make_trgt_sample(
        self,
        sample_name: str,
        region_chrom: str,
        region_start: int,
        region_end: int,
        trid: str,
        sam_records: list[str],
    ) -> Path:
        trgt_prefix = self.root / sample_name
        vcf_text = (
            "##fileformat=VCFv4.2\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n"
            f"{region_chrom}\t{region_start}\t.\tA\t<TR>\t.\tPASS\tTRID={trid};END={region_end}\tGT\t0/1\n"
        )
        self.write_gzip_text(Path(f"{trgt_prefix}.trgt.vcf.gz"), vcf_text)

        sam_text = (
            "@HD\tVN:1.6\tSO:coordinate\n"
            f"@SQ\tSN:{region_chrom}\tLN:1000\n"
            + "".join(sam_records)
        )
        sam_path = Path(f"{trgt_prefix}.trgt.spanning.sam")
        bam_path = Path(f"{trgt_prefix}.trgt.spanning.bam")
        self.write_text(sam_path, sam_text)
        subprocess.run(
            [
                "samtools",
                "view",
                "-b",
                "-o",
                str(bam_path),
                str(sam_path),
            ],
            check=True,
        )
        return trgt_prefix

    def write_sample_table(self, rows: list[dict[str, str]]) -> Path:
        sample_table = self.root / "samples.tsv"
        fieldnames: list[str] = []
        for row in rows:
            for field_name in row:
                if field_name not in fieldnames:
                    fieldnames.append(field_name)
        with sample_table.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
            writer.writeheader()
            writer.writerows(rows)
        return sample_table

    def test_target_region_parses_multiple_formats(self) -> None:
        region_a = TargetRegion.from_region_text("chr1:10-20")
        region_b = TargetRegion.from_region_text("chr1_10_20")
        self.assertEqual(region_a.region, "chr1_10_20")
        self.assertEqual(region_b.igv_locus, "chr1:10-20")

    @unittest.skipUnless(shutil.which("samtools"), "samtools not available")
    def test_run_workflow_creates_outputs(self) -> None:
        sample_a_dir = self.make_medaka_bam_sample(
            sample_name="sampleA",
            ref_headers=[
                ">tr_chr1_10_20_pad_0_30_fwd_hap1_phased-set0_ploidy2\n",
                ">tr_chr1_10_20_pad_0_30_fwd_hap2_phased-set0_ploidy2\n",
            ],
            sam_records=[
                "r1\t0\ttr_chr1_10_20_pad_0_30_fwd_hap1_phased-set0_ploidy2\t1\t60\t4M\t*\t0\t0\tAAAA\t*\n",
                "r2\t0\ttr_chr1_10_20_pad_0_30_fwd_hap2_phased-set0_ploidy2\t1\t60\t6M\t*\t0\t0\tAAAAAA\t*\n",
            ],
        )
        sample_b_dir = self.make_medaka_bam_sample(
            sample_name="sampleB",
            ref_headers=[
                ">tr_HOM_chr1_10_20_pad_0_30_fwd_hap1_phased-set0_ploidy2\n",
            ],
            sam_records=[
                "r3\t0\ttr_HOM_chr1_10_20_pad_0_30_fwd_hap1_phased-set0_ploidy2\t1\t60\t5M\t*\t0\t0\tAAAAA\t*\n",
                "r4\t0\ttr_HOM_chr1_10_20_pad_0_30_fwd_hap1_phased-set0_ploidy2\t1\t60\t7M\t*\t0\t0\tAAAAAAA\t*\n",
            ],
        )

        sample_table = self.write_sample_table(
            [
                {
                    "sample": "sampleA",
                    "medaka_folder": str(sample_a_dir),
                    "software": "medaka",
                    "flank_bp": "0",
                    "batch": "B1",
                },
                {
                    "sample": "sampleB",
                    "medaka_folder": str(sample_b_dir),
                    "software": "medaka",
                    "flank_bp": "0",
                    "batch": "B2",
                },
            ]
        )

        output_dir = self.root / "outputs"
        result = run_workflow(
            sample_table_path=sample_table,
            output_dir=output_dir,
            region_text="chr1:10-20",
            label_columns=["sample", "batch"],
        )

        self.assertEqual(result["read_source"], "bam")
        self.assertEqual(len(result["figure_paths"]), 1)
        self.assertTrue(result["figure_paths"][0].is_file())
        self.assertTrue(result["summary_tsv_path"].is_file())
        self.assertTrue(result["validation_tsv_path"].is_file())
        self.assertTrue(result["reads_tsv_path"].is_file())

        with gzip.open(result["reads_tsv_path"], "rt") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual(len(rows), 4)
        self.assertEqual({row["hap"] for row in rows}, {"hap1", "hap2"})
        self.assertEqual(rows[0]["flank_bp"], "0")

    @unittest.skipUnless(shutil.which("samtools"), "samtools not available")
    def test_run_workflow_uses_inline_metadata_columns(self) -> None:
        sample_a_dir = self.make_medaka_bam_sample(
            sample_name="sampleA",
            ref_headers=[
                ">tr_chr1_10_20_pad_0_30_fwd_hap1_phased-set0_ploidy2\n",
                ">tr_chr1_10_20_pad_0_30_fwd_hap2_phased-set0_ploidy2\n",
            ],
            sam_records=[
                "r1\t0\ttr_chr1_10_20_pad_0_30_fwd_hap1_phased-set0_ploidy2\t1\t60\t4M\t*\t0\t0\tAAAA\t*\n",
                "r2\t0\ttr_chr1_10_20_pad_0_30_fwd_hap2_phased-set0_ploidy2\t1\t60\t6M\t*\t0\t0\tAAAAAA\t*\n",
            ],
        )
        sample_b_dir = self.make_medaka_bam_sample(
            sample_name="sampleB",
            ref_headers=[
                ">tr_chr1_10_20_pad_0_30_fwd_hap1_phased-set0_ploidy2\n",
                ">tr_chr1_10_20_pad_0_30_fwd_hap2_phased-set0_ploidy2\n",
            ],
            sam_records=[
                "r3\t0\ttr_chr1_10_20_pad_0_30_fwd_hap1_phased-set0_ploidy2\t1\t60\t5M\t*\t0\t0\tAAAAA\t*\n",
                "r4\t0\ttr_chr1_10_20_pad_0_30_fwd_hap2_phased-set0_ploidy2\t1\t60\t7M\t*\t0\t0\tAAAAAAA\t*\n",
            ],
        )
        sample_table = self.write_sample_table(
            [
                {
                    "sample": "sampleA",
                    "medaka_folder": str(sample_a_dir),
                    "software": "medaka",
                    "flank_bp": "0",
                    "tissue": "Lung",
                    "technology": "ONT",
                },
                {
                    "sample": "sampleB",
                    "medaka_folder": str(sample_b_dir),
                    "software": "medaka",
                    "flank_bp": "0",
                    "tissue": "Heart",
                    "technology": "PacBio",
                },
            ]
        )

        output_dir = self.root / "metadata_outputs"
        result = run_workflow(
            sample_table_path=sample_table,
            output_dir=output_dir,
            region_text="chr1:10-20",
            label_columns=["tissue", "technology"],
        )

        self.assertEqual(result["read_source"], "bam")
        with open(result["summary_tsv_path"]) as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual(rows[0]["tissue"], "Lung")
        self.assertEqual(rows[0]["technology"], "ONT")

    @unittest.skipUnless(shutil.which("samtools"), "samtools not available")
    def test_flank_bp_is_subtracted_from_read_length(self) -> None:
        sample_a_dir = self.make_medaka_bam_sample(
            sample_name="sampleA",
            ref_headers=[">tr_chr1_10_20_pad_0_30_fwd_hap1_phased-set0_ploidy2\n"],
            sam_records=[
                "r1\t0\ttr_chr1_10_20_pad_0_30_fwd_hap1_phased-set0_ploidy2\t1\t60\t24M\t*\t0\t0\t"
                + ("A" * 24)
                + "\t*\n"
            ],
        )
        sample_table = self.write_sample_table(
            [
                {
                    "sample": "sampleA",
                    "medaka_folder": str(sample_a_dir),
                    "software": "medaka",
                    "flank_bp": "5",
                }
            ]
        )

        output_dir = self.root / "padding_outputs"
        result = run_workflow(
            sample_table_path=sample_table,
            output_dir=output_dir,
            region_text="chr1:10-20",
        )

        with gzip.open(result["reads_tsv_path"], "rt") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0]["raw_read_length"], "24")
        self.assertEqual(rows[0]["read_length"], "14")
        self.assertEqual(rows[0]["flank_bp"], "5")

    def test_run_workflow_requires_software_and_flank_bp_columns(self) -> None:
        sample_table = self.write_sample_table(
            [
                {
                    "sample": "sampleA",
                    "medaka_folder": "/tmp/sampleA.medaka",
                },
            ]
        )

        with self.assertRaisesRegex(ValueError, "Missing: .*flank_bp.*software|Missing: .*software.*flank_bp"):
            run_workflow(
                sample_table_path=sample_table,
                output_dir=self.root / "invalid_outputs",
                region_text="chr1:10-20",
            )

    @unittest.skipUnless(shutil.which("samtools"), "samtools not available")
    def test_run_workflow_supports_mixed_medaka_and_trgt_inputs(self) -> None:
        medaka_dir = self.make_medaka_bam_sample(
            sample_name="sampleA",
            ref_headers=[
                ">tr_chr1_10_20_pad_0_30_fwd_hap1_phased-set0_ploidy2\n",
                ">tr_chr1_10_20_pad_0_30_fwd_hap2_phased-set0_ploidy2\n",
            ],
            sam_records=[
                "r1\t0\ttr_chr1_10_20_pad_0_30_fwd_hap1_phased-set0_ploidy2\t1\t60\t8M\t*\t0\t0\tAAAAAAAA\t*\n",
            ],
        )
        trgt_prefix = self.make_trgt_sample(
            sample_name="sampleT",
            region_chrom="chr1",
            region_start=10,
            region_end=20,
            trid="12",
            sam_records=[
                "tr1\t0\tchr1\t10\t60\t10M\t*\t0\t0\tAAAAAAAAAA\t*\tTR:Z:12\tHP:i:1\n",
                "tr2\t0\tchr1\t10\t60\t12M\t*\t0\t0\tAAAAAAAAAAAA\t*\tTR:Z:12\tHP:i:2\n",
                "off_target\t0\tchr1\t10\t60\t14M\t*\t0\t0\tAAAAAAAAAAAAAA\t*\tTR:Z:99\tHP:i:1\n",
            ],
        )

        sample_table = self.write_sample_table(
            [
                {
                    "sample": "sampleA",
                    "medaka_folder": str(medaka_dir),
                    "software": "medaka",
                    "flank_bp": "1",
                },
                {
                    "sample": "sampleT",
                    "medaka_folder": str(trgt_prefix),
                    "software": "trgt",
                    "flank_bp": "2",
                },
            ]
        )

        output_dir = self.root / "mixed_outputs"
        result = run_workflow(
            sample_table_path=sample_table,
            output_dir=output_dir,
            region_text="chr1:10-20",
            label_columns=["sample", "software"],
        )

        self.assertEqual(result["read_source"], "bam")
        with gzip.open(result["reads_tsv_path"], "rt") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))

        self.assertEqual(len(rows), 3)
        rows_by_sample = {row["sample"]: [] for row in rows}
        for row in rows:
            rows_by_sample[row["sample"]].append(row)

        medaka_row = rows_by_sample["sampleA"][0]
        self.assertEqual(medaka_row["software"], "medaka")
        self.assertEqual(medaka_row["read_length"], "6")
        self.assertEqual(medaka_row["flank_bp"], "1")

        trgt_lengths = sorted(int(row["read_length"]) for row in rows_by_sample["sampleT"])
        self.assertEqual(trgt_lengths, [6, 8])
        self.assertEqual({row["software"] for row in rows_by_sample["sampleT"]}, {"trgt"})
        self.assertEqual({row["hap"] for row in rows_by_sample["sampleT"]}, {"hap1", "hap2"})
        self.assertEqual({row["flank_bp"] for row in rows_by_sample["sampleT"]}, {"2"})

    @unittest.skipUnless(shutil.which("samtools"), "samtools not available")
    def test_run_workflow_supports_multi_region_trgt_bed_input(self) -> None:
        trgt_prefix = self.root / "sampleT"
        vcf_text = (
            "##fileformat=VCFv4.2\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n"
            "chr1\t10\t.\tA\t<TR>\t.\tPASS\tTRID=12;END=20\tGT\t0/1\n"
            "chr1\t30\t.\tA\t<TR>\t.\tPASS\tTRID=34;END=40\tGT\t0/1\n"
        )
        self.write_gzip_text(Path(f"{trgt_prefix}.trgt.vcf.gz"), vcf_text)

        sam_text = (
            "@HD\tVN:1.6\tSO:coordinate\n"
            "@SQ\tSN:chr1\tLN:1000\n"
            "tr1\t0\tchr1\t10\t60\t10M\t*\t0\t0\tAAAAAAAAAA\t*\tTR:Z:12\tHP:i:1\n"
            "tr2\t0\tchr1\t30\t60\t12M\t*\t0\t0\tAAAAAAAAAAAA\t*\tTR:Z:34\tHP:i:2\n"
            "off_target\t0\tchr1\t30\t60\t8M\t*\t0\t0\tAAAAAAAA\t*\tTR:Z:99\tHP:i:1\n"
        )
        sam_path = Path(f"{trgt_prefix}.trgt.spanning.sam")
        bam_path = Path(f"{trgt_prefix}.trgt.spanning.bam")
        self.write_text(sam_path, sam_text)
        subprocess.run(
            [
                "samtools",
                "view",
                "-b",
                "-o",
                str(bam_path),
                str(sam_path),
            ],
            check=True,
        )
        sorted_bam_path = Path(f"{trgt_prefix}.trgt.spanning.sorted.bam")
        subprocess.run(
            [
                "samtools",
                "sort",
                "-o",
                str(sorted_bam_path),
                str(bam_path),
            ],
            check=True,
        )
        subprocess.run(["samtools", "index", str(sorted_bam_path)], check=True)

        sample_table = self.write_sample_table(
            [
                {
                    "sample": "sampleT",
                    "medaka_folder": str(trgt_prefix),
                    "software": "trgt",
                    "flank_bp": "2",
                },
            ]
        )
        bed_path = self.root / "targets.bed"
        self.write_text(bed_path, "chr1\t10\t20\nchr1\t30\t40\n")

        output_dir = self.root / "multi_region_trgt_outputs"
        result = run_workflow(
            sample_table_path=sample_table,
            output_dir=output_dir,
            bed_path=bed_path,
            label_columns=["sample", "software"],
        )

        with gzip.open(result["reads_tsv_path"], "rt") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))

        self.assertEqual(len(rows), 2)
        self.assertEqual({row["region"] for row in rows}, {"chr1_10_20", "chr1_30_40"})
        self.assertEqual({row["hap"] for row in rows}, {"hap1", "hap2"})
        temp_dir = output_dir / ".tmp"
        self.assertTrue(temp_dir.is_dir())
        self.assertEqual(list(temp_dir.iterdir()), [])


if __name__ == "__main__":
    unittest.main()
