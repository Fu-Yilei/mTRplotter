# mTRplotter

`mTRplotter` turns the read-length plotting logic from the notebook into a small Python CLI.

It expects a tab-delimited sample table with at least:

- `sample`
- `medaka_folder`
- `software` (`medaka` or `trgt`)
- `flank_bp`
- any number of extra metadata columns

For each requested locus, it:

- validates medaka rows against each sample's `ref_chunks.fasta`
- validates TRGT rows against each sample's `.trgt.vcf.gz`
- extracts medaka reads from `trimmed_reads_to_poa.bam` through the BAM index
- extracts TRGT reads from `.trgt.spanning.bam`
- extracts read lengths for the requested locus or BED subset
- writes per-read and summary tables
- produces one figure per locus with samples on the x-axis and trimmed read length on the y-axis

## Usage

Single locus:

```bash
python -m mtrplotter \
  --sample-table examples/smht001_ont_pacbio_samples.tsv \
  --region chr1:99682978-99683197 \
  --output-dir demo_outputs/smht001_ont_pacbio_chr1_99682978_99683197_inline_metadata_jobs10 \
  --label-columns tissue_name tech core technology \
  --jobs 10
```

BED input:

```bash
python -m mtrplotter \
  --sample-table samples.tsv \
  --bed regions.bed \
  --output-dir results \
  --catalog-bed catalog.bed
```

If `--catalog-bed` is not provided, loci are still checked against each sample's `ref_chunks.fasta`, which reflects the regions actually used by the medaka run.

Metadata should already be present in the sample table. There is no separate manifest join step at runtime.

The bundled example table keeps inline columns such as `tissue`, `tissue_name`, `tech`, `core`, and `technology` so plotting labels can be chosen directly from the sample table.

For medaka rows, the plotted length is `raw_length - 2 * flank_bp`, where `flank_bp` comes from the sample table.

For TRGT rows, the software looks for the `FL` BAM tag and uses `raw_length - left_flank - right_flank`. If the `FL` tag is absent, it falls back to the sample-table `flank_bp`.

For TRGT rows, use the `medaka_folder` column as a generic input path. The most reliable form is a per-sample TRGT prefix such as:

```tsv
sample	medaka_folder	software	flank_bp
SMHT001-3A-XX-M42-B001-broad	/path/to/trgt/SMHT001-3A-XX-M42-B001-broad	trgt	50
```

`mTRplotter` resolves the sibling files:

- `.trgt.vcf.gz`
- `.trgt.spanning.bam`
