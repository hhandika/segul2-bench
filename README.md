# segul2-bench

Benchmark SEGUL compared to other software.

## Setup

Python is managed with [uv](https://docs.astral.sh/uv/).

```bash
uv sync
```

Tools: SEGUL is built from the working tree at `~/Codes/Rust/segul`
(`cargo build --release`, done by the harness). fastp is read from
`~/Programs/fastp` and FastQC from `apps/FastQC/`. Paths are set in
`bench/*.toml`. Datasets live in `data/`.

## Read summary: SEGUL vs fastp and FastQC

Compares `segul read summary` (minimal, default, complete) on the gzipped
FASTQ files in `data/` with:

- fastp in report-only mode (`-A -G -Q -L`, no output reads).
- FastQC running only the modules SEGUL complete also computes (per-tile
  quality and k-mer content ignored), so both compute identical features.

Wall time and peak memory come from `/usr/bin/time -v`, with matching thread
counts.

```bash
uv run segul2-bench read-summary list    # show the commands
uv run segul2-bench read-summary run     # run (results/read-summary/latest.csv)
uv run jupyter nbconvert --to notebook --execute --inplace notebooks/read_summary.ipynb
```

Figures are saved to `figures/read-summary/`. See `bench/read-summary.toml`
for settings and for what each tool computes in each mode.
