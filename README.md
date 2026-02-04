# asm-smk

A production-ready pipeline for generating high-quality, phased diploid genome assemblies from PacBio HiFi reads. Supports trio-based and Hi-C phasing, optional ultra-long ONT reads, and produces reference-oriented assemblies with standardized headers ready for pangenome analysis.

## Install

Please start by installing [pixi](https://pixi.sh/latest/) which handles the environment of this Snakemake workflow.

You can then install the `pixi` environment by cloning this repository and running:

```bash
pixi install
```

## Test case

Before running the workflow please run the test case to make sure everything is working as expected.

```bash
pixi run test
```

## Usage

`pixi` handles the execution of the Snakemake workflows:

```bash
pixi run snakemake --configfile ...
```

And if you want to run this Snakemake from another directory you can do so with:

```bash
pixi run --manifest-path /path/to/snakemake/pixi.toml snakemake ...
```

where you update `/path/to/snakemake/pixi.toml` to the path of the `pixi.toml` you cloned.

And in place of `...` use all the normal Snakemake arguments for your workflow.

## Configuration

Your configuration file (e.g. `config.yaml`) should have a manifest entry that points to the inputs to be assembled. For example:

```yaml
manifest: test/test.tbl
```

### Manifest Format

The manifest file should be a whitespace-separated file with the following columns:

| Column      | Required | Description                                      |
| ----------- | -------- | ------------------------------------------------ |
| `sample`    | Yes      | Unique sample identifier                         |
| `hifi`      | Yes      | Path(s) to HiFi reads                            |
| `paternal`  | No       | Path(s) to paternal short reads for trio phasing |
| `maternal`  | No       | Path(s) to maternal short reads for trio phasing |
| `hic_r1`    | No       | Path(s) to Hi-C R1 reads for Hi-C phasing        |
| `hic_r2`    | No       | Path(s) to Hi-C R2 reads for Hi-C phasing        |
| `ultralong` | No       | Path(s) to ultra-long ONT reads                  |

The inputs need not be in fastq format: bam, sam, cram, and fasta are also supported. The workflow will automatically detect the file type and run the appropriate tools.

Use `NA` for columns where you don't have data.

### Assembly Modes

The workflow automatically determines the assembly mode based on provided inputs:

| Mode  | Phasing | Required Columns               | Optional    |
| ----- | ------- | ------------------------------ | ----------- |
| `bp`  | None    | `hifi`                         | `ultralong` |
| `dip` | Trio    | `hifi`, `paternal`, `maternal` | `ultralong` |
| `hic` | Hi-C    | `hifi`, `hic_r1`, `hic_r2`     | `ultralong` |

**Important:** Trio phasing and Hi-C phasing are **mutually exclusive**. Do not provide both parental data and Hi-C data for the same sample.

### Example Manifest

```
sample       hifi               paternal      maternal      hic_r1          hic_r2          ultralong
HG002_bp     /data/hifi.fq.gz   NA            NA            NA              NA              NA
HG002_trio   /data/hifi.fq.gz   /pat.fq.gz    /mat.fq.gz    NA              NA              NA
HG002_hic    /data/hifi.fq.gz   NA            NA            /hic_R1.fq.gz   /hic_R2.fq.gz   NA
HG002_hic_ul /data/hifi.fq.gz   NA            NA            /hic_R1.fq.gz   /hic_R2.fq.gz   /ul.fq.gz
HG002_dip_ul /data/hifi.fq.gz   /pat.fq.gz    /mat.fq.gz    NA              NA              /ul.fq.gz
```

### Multiple Input Files

All read columns support comma-separated lists for multiple input files:

```
sample  hifi                   hic_r1                  hic_r2                  ultralong
multi   /f1.fq.gz,/f2.fq.gz    /r1a.fq.gz,/r1b.fq.gz   /r2a.fq.gz,/r2b.fq.gz   /ul1.fq,/ul2.fq
```

## Output

The workflow produces oriented assemblies and per-reference alignments:

```
results/
├── assemblies/
│   ├── {sample}.{mode}.hap1.fa.gz      # Haplotype 1 assembly (oriented to primary ref)
│   ├── {sample}.{mode}.hap1.fa.gz.fai  # FASTA index
│   ├── {sample}.{mode}.hap2.fa.gz      # Haplotype 2 assembly
│   ├── {sample}.{mode}.hap2.fa.gz.fai
│   └── {sample}.{mode}.dip.fa.gz       # Combined diploid assembly
├── alignments/
│   ├── {reference}/
│   │   ├── {sample}.{mode}.hap1.bam        # Haplotype 1 aligned to reference
│   │   ├── {sample}.{mode}.hap1.bam.csi    # BAM index
│   │   ├── {sample}.{mode}.hap1.paf        # PAF format alignment
│   │   ├── {sample}.{mode}.hap2.bam
│   │   ├── {sample}.{mode}.hap2.paf
│   │   └── {sample}.{mode}.dip.bam         # Combined diploid alignment
│   └── {reference2}/
│       └── ...
└── temp/                                   # Intermediate files (auto-deleted)
```

Assemblies are oriented once to the primary reference (T2T-CHM13v2.0 by default), while alignments are produced for all configured references.

Where `{mode}` is one of: `bp` (unphased), `dip` (trio-phased), or `hic` (Hi-C phased).

### PanSN-spec Headers

Assembly FASTA headers follow the [PanSN-spec](https://github.com/pangenome/PanSN-spec) naming convention for compatibility with pangenome tools:

```
>{sample}#{haplotype}#{contig} assembler=hifiasm version=X.X.X phasing={method} ultralong={yes|no}
```

Example:

```
>HG002#1#ptg000001l assembler=hifiasm version=0.25.0 phasing=trio ultralong=yes
>HG002#2#ptg000001l assembler=hifiasm version=0.25.0 phasing=trio ultralong=yes
```

Header metadata fields:

- `assembler` - Assembler used (e.g., `hifiasm`)
- `version` - Assembler version
- `phasing` - Phasing method: `trio`, `hic`, or `none`
- `ultralong` - Whether ultra-long ONT reads were used: `yes` or `no`
- `orient_ref` - Reference used for orientation (only in reference-oriented assemblies)
- `chrom` - Primary chromosome the contig aligns to (only in reference-oriented assemblies)
- `flipped` - Whether the contig was reverse complemented: `yes` or `no` (only in reference-oriented assemblies)

### Reference-Oriented Assemblies

Assemblies are oriented to the **primary reference** (first in config, T2T-CHM13v2.0 by default). Contigs are flipped to match reference strand orientation and annotated with chromosome assignments.

The orientation pipeline:
1. **Fast alignment with minimap2** - Quick whole-genome alignment to determine contig orientation
2. **Flip contigs** - Reverse complement contigs where majority of aligned bases are on reverse strand
3. **Full alignment with minimap2** - Detailed alignment to all configured references

Reference-oriented assemblies have additional header annotations:
- `orient_ref` - Reference used for orientation (e.g., `T2T-CHM13v2.0`)
- `chrom` - Primary chromosome the contig aligns to (e.g., `chr1`)
- `flipped` - Whether the contig was reverse complemented (`yes` or `no`)

Example oriented header:

```
>HG002#1#ptg000001l assembler=hifiasm version=0.25.0 phasing=trio ultralong=yes orient_ref=T2T-CHM13v2.0 chrom=chr1 flipped=no
>HG002#1#ptg000015l assembler=hifiasm version=0.25.0 phasing=trio ultralong=yes orient_ref=T2T-CHM13v2.0 chrom=chr7 flipped=yes
```

To configure references, add them to your config file:

```yaml
manifest: samples.tbl
references:
  T2T-CHM13v2.0: /path/to/chm13v2.fa   # Primary (used for orientation)
  GRCh38: /path/to/grch38.fa           # Alignments only
```

When no reference is provided, default references (T2T-CHM13v2.0 and GRCh38) are downloaded automatically.

### Submitting to the Hyak HPC via Slurm

```bash
pixi run snakemake --configfile /path/to/your/config.yaml --profile profiles/slurm-executor
```
