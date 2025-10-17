# rsv_pipeline

Genomic surveillance workflow for Respiratory Syncytial Virus (RSV) built on Snakemake and packaged for installation via the [`snk`](https://snakemake.github.io/snk) workflow manager.

## Prerequisites

- [Mamba](https://mamba.readthedocs.io/en/latest/installation.html) or Conda for environment management.
- The `snk` CLI (installed below).
- Access to paired-end FASTQ reads produced with ARTIC-style schemes.

## Installation with `snk`

```bash
# create and activate a minimal environment that ships with snk
mamba create -n snk -c bioconda snk
mamba activate snk

# install the workflow (ships as the CLI command `rsvp`)
snk install charlesfoster/rsv_pipeline --name rsvp --force --dependency pandas
```

The command above places the workflow and its Snakemake environments under the active `snk` environment. Update it later with `snk install … --force`.

## Running the workflow

```bash
rsvp run \
  --reads-dir /path/to/reads \
  --outdir results/rsv_run_01 \
  --cores 8 \
  --suffix1 _L001_R1_001.fastq.gz \
  --suffix2 _L001_R2_001.fastq.gz
```

Key convenience flags include:

- `--reads-dir` / `--outdir`: required input and output roots.
- `--remove-sample-index` / `--no-remove-sample-index`: control how sample names are parsed from FASTQ filenames.
- `--keep-trimmed-reads`: persist the intermediate FASTQ files if needed.
- `--wastewater-mode`: enable the Freyja demixing and plotting branch.
- `--freyja-sampling-times`, `--freyja-plot-interval`, `--freyja-custom-lineages`: optional metadata for Freyja plots.
- `--min-depth`, `--consensus-threshold`, `--snp-threshold`, `--indel-threshold`, `--mixed-threshold`: adjust variant calling heuristics without editing `config.yaml`.

The `rsvp run` wrapper forwards any unknown flags to Snakemake, so standard options such as `--quiet`, `--dry`, `--dag report.pdf`, and `--profile myprofile` behave exactly as they do in bare Snakemake. 

To inspect every available option, run:

```bash
rsvp run --help
```

### First run
When _first_ running the pipeline, it's necessary to copy the data and scripts from the `snk` installation like so:

```
rsvp run \
  -r data \
  -r scripts \
  --reads-dir /path/to/reads/folder \
  --outdir results/test \
  --keep-resources 
```

On subsequent runs, the `-r data`, `-r scripts` and `--keep-resources` parts can be omitted if:
- you are running the pipeline from the same directory each time
- you do not need to update the `data` and `scripts` folder, e.g. if the pipeline has changed

### Rerunning specific steps

Because `rsvp run` simply wraps Snakemake, you can bring forward any rerun behaviour by passing through native Snakemake switches. For example, to recompute Freyja demixing outputs:

```bash
rsvp run \
  --reads-dir /path/to/reads/folder \
  --outdir results/rsv_run_01 \
  --wastewater-mode \
  --forcerun freyja_demix
```

### Configuration file (optional)

If you prefer working with a YAML configuration file, create one based on `config.yaml` and point to it with `--config my_config.yaml`. CLI flags always take precedence over values defined in the file.

## Troubleshooting & tips

- Ensure `suffix1`/`suffix2` match the read naming pattern; otherwise no samples will be detected.
- Use `--dry` to preview the DAG, or `--dag dag.pdf` to visualise the full graph.
- Add `--keep-resources` or `--keep-snakemake` if you need to inspect intermediate databases or Snakemake bookkeeping after a run.
- For reproducible deployments (HPC, cloud), create a Snakemake profile and supply it with `--profile`.
