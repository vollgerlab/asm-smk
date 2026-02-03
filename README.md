# asm-smk

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

| Column | Required | Description |
|--------|----------|-------------|
| `sample` | Yes | Unique sample identifier |
| `hifi` | Yes | Path(s) to HiFi reads |
| `paternal` | No | Path(s) to paternal short reads for trio phasing |
| `maternal` | No | Path(s) to maternal short reads for trio phasing |
| `hic_r1` | No | Path(s) to Hi-C R1 reads for Hi-C phasing |
| `hic_r2` | No | Path(s) to Hi-C R2 reads for Hi-C phasing |
| `ultralong` | No | Path(s) to ultra-long ONT reads |

The inputs need not be in fastq format: bam, sam, cram, and fasta are also supported. The workflow will automatically detect the file type and run the appropriate tools.

Use `NA` for columns where you don't have data.

### Assembly Modes

The workflow automatically determines the assembly mode based on provided inputs:

| Mode | Phasing | Required Columns | Optional |
|------|---------|------------------|----------|
| `bp` | None | `hifi` | `ultralong` |
| `dip` | Trio | `hifi`, `paternal`, `maternal` | `ultralong` |
| `hic` | Hi-C | `hifi`, `hic_r1`, `hic_r2` | `ultralong` |

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

### Submitting to the Hyak HPC via Slurm

```bash
pixi run snakemake --configfile /path/to/your/config.yaml --profile profiles/slurm-executor
```
