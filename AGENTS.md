# Agent instructions

segul2-bench benchmarks SEGUL against other phylogenomics and read QC
software (fastp, FastQC, ...) on real datasets. This file tells coding agents
how the repository is organized and how to work in it.

## Hard rules

- Never commit. Do not run `git commit`, `git push`, or anything that creates
  commits or tags. Leave all changes uncommitted for the maintainer to review
  and commit.
- Do not publish packages (`uv publish`), create GitHub releases, or push
  tags.
- Do not delete or overwrite files you did not create in this task. This
  includes datasets in `data/`, tools in `apps/`, earlier results in
  `results/`, and anything outside the repo (the SEGUL repo, test datasets).
- Python is always managed with `uv`. Add dependencies with `uv add` (or
  `uv add --dev`), run code with `uv run ...`, and install with `uv sync`.
  Never use `pip`, `python -m pip`, a bare `python`, or a manually created
  virtualenv.

## Repository layout

- `pyproject.toml`, `uv.lock`, `.python-version`: the uv project
  (`segul2-bench`, Python 3.12). The `segul2-bench` command is defined under
  `[project.scripts]`.
- `src/segul2_bench/`: the benchmark harness.
  - `cli.py`: command-line entry point (`segul2-bench <benchmark> <action>`).
  - `config.py`: loads a benchmark TOML config into dataclasses.
  - `timing.py`: runs one command under `/usr/bin/time -v` and parses wall
    time, CPU time, and peak memory (max RSS).
  - `read_summary.py`: the read summary benchmark (SEGUL `read summary` vs
    fastp and FastQC).
- `bench/`: benchmark configs, one TOML per benchmark
  (`bench/read-summary.toml`). Tool paths, thread counts, repeats, and
  datasets are set here, not in code.
- `notebooks/`: Jupyter notebooks that load results and plot them
  (`read_summary.ipynb`). Figures are written to `figures/<benchmark>/`.
- `data/`: input datasets (ignored by git, large). Benchmarks always use the
  data in this folder, not datasets elsewhere on the machine.
  - `data/raw-reads/`: 3 paired raw gzipped FASTQ samples (`*_READ1.fq.gz`,
    `*_READ2.fq.gz`, about 0.9 GB).
  - `data/alignments/Bunomys_*/`: 5 paired gzipped FASTQ samples of reads
    already cleaned by fastp (`*_R1.fq.gz`, `*_R2.fq.gz`, about 1.9 GB), with
    the fastp cleaning reports in `reports/`. This is the SEGUL `read-tests`
    set. The folder may move as `data/` is reorganized; update
    `bench/read-summary.toml` if it does.
  - `data/alignments/<study>/`: published alignment datasets (NEXUS loci and
    concatenated alignments) for alignment benchmarks.
  - `data/cleaned-reads/`: reserved for cleaned reads.
- `apps/`: third-party tools installed locally (ignored by git). FastQC is
  `apps/FastQC/fastqc` (v0.13, runs on the system Java). fastp is the static
  binary at `~/Programs/fastp`. Point the config to the binaries rather than
  assuming they are on `PATH`.
- `results/`: raw benchmark output (ignored by git).
  `results/<benchmark>/<date>-<HHMM>-<host>.csv` holds one row per run, with a
  matching `.json` holding machine and tool metadata, and `latest.csv` /
  `latest.json` copies of the newest run. `results/<benchmark>/work/` is
  scratch output from the tools and can be removed by the harness.
- `figures/`: plots produced by the notebooks.

## Running benchmarks

```bash
uv sync
uv run segul2-bench read-summary run                  # full run, config defaults
uv run segul2-bench read-summary run --repeats 1 --dataset raw-reads   # quick check
uv run segul2-bench read-summary list                 # show jobs without running
```

- SEGUL is the working-tree release build from `segul_repo` in the config
  (`~/Codes/Rust/segul`). The harness runs `cargo build --release` there and
  copies the binary to `.cache/bin/` so a rebuild during a run does not change
  the binary being measured. Pass `--no-build` to use the existing build.
- Full runs take a long time (about 30 minutes; fastp and FastQC are the slow
  part). Run them in the
  background and redirect output to a log in `results/`.
- Each command runs once as an unrecorded warmup to fill the page cache, then
  `repeats` times. Tools are interleaved within a repeat so drift affects
  both equally.
- Thread counts match between tools: SEGUL gets `RAYON_NUM_THREADS`, fastp
  gets `-w`, FastQC gets `--threads`. FastQC uses at most 4 threads per file,
  so single-file FastQC runs are paired with a SEGUL complete run at the same
  thread count. Keep this when adding tools.
- Keep comparisons fair and documented: run other tools in report-only mode
  where possible, and record in the config (`computes` fields) what each tool
  computes in each mode. The notebook prints these notes next to the plots.
- Compare only identical features where the tool allows it. FastQC modules
  SEGUL does not compute (per-tile quality, k-mer content) are listed under
  `[tools.fastqc] ignore` and switched off through a generated limits file
  (`results/read-summary/work/fastqc-limits.txt`).

## Notebooks

```bash
uv run jupyter lab                                     # interactive
uv run jupyter nbconvert --to notebook --execute --inplace notebooks/read_summary.ipynb
```

- Notebooks read `results/<benchmark>/latest.csv`, never re-run tools.
- Plotting uses pandas and matplotlib. Save figures to `figures/<benchmark>/`
  as PNG.
- Colors are ColorBrewer Dark2 (`matplotlib.colormaps["Dark2"]`), one fixed
  color per tool: SEGUL slot 1 (teal), fastp slot 2 (orange), FastQC slot 3
  (purple). New tools take the next slot. Modes are told apart by marker and
  line style, not color.
- Commit notebooks executed, so plots show on GitHub.

## Coding style

- Python 3.12 with type hints, short docstrings, and `pathlib.Path` for paths.
- Harness code (`src/`) uses only the standard library. Analysis
  dependencies (pandas, matplotlib, jupyter) are for notebooks.
- Format with `uv run ruff format` and check with `uv run ruff check` when
  ruff is installed as a dev dependency.
- Tool invocations are built as argument lists, never shell strings.
- New benchmarks get a TOML in `bench/`, a module in `src/segul2_bench/`, a
  subcommand in `cli.py`, a notebook in `notebooks/`, and a section in this
  file and in `README.md`.
