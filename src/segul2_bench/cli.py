"""Command-line entry point: `uv run segul2-bench <benchmark> <action>`."""

import argparse
import shutil
import socket
from datetime import datetime
from pathlib import Path

from segul2_bench.config import REPO_ROOT, load_read_summary
from segul2_bench.read_summary import (
    ReadSummaryBench,
    build_segul,
    freeze_segul,
    write_metadata,
)


def main() -> None:
    parser = argparse.ArgumentParser(
        prog="segul2-bench", description="Benchmark SEGUL against other software"
    )
    benches = parser.add_subparsers(dest="bench", required=True)

    reads = benches.add_parser(
        "read-summary", help="SEGUL read summary vs fastp and FastQC"
    )
    reads.add_argument("action", choices=["run", "list"])
    reads.add_argument(
        "--config", type=Path, default=REPO_ROOT / "bench" / "read-summary.toml"
    )
    reads.add_argument("--dataset", action="append", help="Dataset id (repeatable)")
    reads.add_argument(
        "--scope", action="append", choices=["file", "dataset"], help="Repeatable"
    )
    reads.add_argument(
        "--tool",
        action="append",
        choices=["segul", "fastp", "fastqc"],
        help="Repeatable",
    )
    reads.add_argument("--repeats", type=int, help="Override config repeats")
    reads.add_argument("--warmup", type=int, help="Override config warmup runs")
    reads.add_argument("--threads", type=int, help="Override config threads")
    reads.add_argument(
        "--no-build", action="store_true", help="Use the existing segul release build"
    )
    reads.add_argument(
        "--output", type=Path, default=REPO_ROOT / "results" / "read-summary"
    )

    args = parser.parse_args()
    if args.bench == "read-summary":
        run_read_summary(args)


def run_read_summary(args: argparse.Namespace) -> None:
    config = load_read_summary(args.config)
    if args.threads:
        config.threads = args.threads
    repeats = args.repeats or config.repeats
    warmup = config.warmup if args.warmup is None else args.warmup
    datasets = [d for d in config.datasets if not args.dataset or d.id in args.dataset]
    scopes = args.scope or ["file", "dataset"]
    tools = args.tool or ["segul", "fastp", "fastqc"]

    segul_bin = config.segul_repo / "target" / "release" / "segul"
    bench = ReadSummaryBench(config, args.output, segul_bin)
    groups = bench.jobs(datasets, scopes, tools)

    if args.action == "list":
        for group in groups:
            for job in group:
                print(
                    f"{job.dataset:10} {job.scope:7} {job.input:45} "
                    f"{job.input_mb:8.1f} MB  {job.label} ({job.threads} threads)"
                )
                for cmd in job.commands:
                    print("    " + " ".join(cmd.argv))
        return

    if "segul" in tools:
        bench.segul_bin = (
            freeze_segul(config.segul_repo)
            if args.no_build
            else build_segul(config.segul_repo)
        )
        groups = bench.jobs(datasets, scopes, tools)

    args.output.mkdir(parents=True, exist_ok=True)
    stem = f"{datetime.now():%Y-%m-%d-%H%M}-{socket.gethostname()}"
    csv_path = args.output / f"{stem}.csv"
    meta_path = args.output / f"{stem}.json"
    write_metadata(meta_path, config, bench.segul_bin, groups, repeats, warmup)
    bench.run(groups, repeats, warmup, csv_path)
    shutil.copy2(csv_path, args.output / "latest.csv")
    shutil.copy2(meta_path, args.output / "latest.json")
    print(f"Results: {csv_path}")
