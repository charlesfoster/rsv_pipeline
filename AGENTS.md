# Repository Guidelines

## Project Structure & Module Organization
- `Snakefile` orchestrates the Snakemake workflow; keep new rules grouped by pipeline stage.
- `config.yaml` holds runtime inputs, thresholds, and suffix parsing options; review before each run.
- `scripts/` contains helper Python and shell utilities consumed by rules; keep filenames snake_case.
- `envs/` defines per-rule Conda env specs invoked with `--use-conda`.
- `data/` ships reference genomes, primers, and the `rsv_database.syldb`; treat as read-only.
- `reads/` and `results/` are default input/output staging areas; avoid committing generated FASTQs or large artefacts.

## Build, Test, and Development Commands
- `mamba env create -f environment.yml` provisions the base development environment.
- `mamba activate rsv_pipeline` activates that environment for local work.
- `snakemake --use-conda --cores 8` runs the full pipeline with per-rule environments; tune `--cores` to available CPUs.
- `snakemake -n --summary` performs a dry run and lists planned outputs; run after any rule or config change.
- `snakemake --unlock` is available if an interrupted run leaves the workspace locked; use cautiously.

## Coding Style & Naming Conventions
- Python scripts use 4-space indentation, clear docstrings, and snake_case functions (`scripts/process_rsv.py` is the model).
- Snakemake rules stay lowercase with underscores (`rule trim_reads`), and helper functions live near the top of `Snakefile`.
- Generated outputs should follow `<sample>/<sample>.<suffix>` under `results/`; reuse that pattern when adding targets.

## Testing Guidelines
- Start with `snakemake -n` to verify DAG integrity, then execute a focused target (e.g. `snakemake results/sample1/sample1.consensus.fasta`) before full runs.
- Inspect `results/overall_samples_qc_summary.csv` and merged Nextclade TSV outputs whenever rules change to confirm column order and counts.
- When adding scripts, include lightweight sanity checks (argument parsing, I/O) that can be run with `python scripts/<name>.py --help`.

## Commit & Pull Request Guidelines
- Git history uses short, imperative subjects (`Added envs for --use-conda`); keep summaries ≤72 chars and expand detail in the body.
- One logical change per commit; reference tickets or samples affected in the body.
- Pull requests should outline the motivation, list validation commands (e.g. `snakemake --use-conda --cores 8`), and attach relevant log snippets or result deltas.

## Configuration Notes
- Keep `config.yaml` synchronized with documentation; note any required `suffix1/suffix2` changes for new sequencing runs.
- Store credentials and raw reads outside the repo; point to them via `reads_dir` in the config rather than hard-coding paths.
