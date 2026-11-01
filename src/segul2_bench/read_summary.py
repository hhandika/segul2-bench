"""Benchmark SEGUL `read summary` against fastp and FastQC.

fastp runs in report-only mode. FastQC runs only the modules SEGUL complete
mode also computes (per-tile quality and k-mer content are ignored through a
limits file).

Two scopes are measured for each dataset:

- `file`: each FASTQ file on its own. SEGUL gets one input, fastp runs in
  single-end mode. Shows how time and memory scale with file size. FastQC
  uses at most 4 threads per file, so SEGUL complete also runs with that many
  threads here to match it.
- `dataset`: all files at once. SEGUL and FastQC summarize every file in one
  call. fastp runs once per R1/R2 pair in paired-end mode, one pair after
  another, which is how fastp is normally used. Time and CPU are summed over
  pairs and peak memory is the largest pair.
"""

import csv
import json
import os
import platform
import shutil
import subprocess
import zipfile
from dataclasses import asdict, dataclass, field
from datetime import datetime
from pathlib import Path

from segul2_bench.config import REPO_ROOT, Dataset, ReadSummaryConfig
from segul2_bench.timing import Measurement, measure

SEGUL_MODES = ("minimal", "default", "complete")
FASTP_MODES = ("report", "complete")
FASTP_BASE_ARGS = ["-A", "-G", "-Q", "-L"]
FASTP_MODE_ARGS = {
    "report": ["--dont_eval_duplication"],
    "complete": ["-p"],
}
# FastQC processes one file with at most this many threads.
FASTQC_THREADS_PER_FILE = 4

FIELDS = [
    "dataset",
    "scope",
    "input",
    "n_files",
    "input_mb",
    "tool",
    "mode",
    "label",
    "threads",
    "repeat",
    "wall_s",
    "user_s",
    "sys_s",
    "max_rss_mb",
    "reads",
    "exit_code",
]


@dataclass
class Command:
    """One process to run, with its own working directory."""

    argv: list[str]
    workdir: Path
    env: dict[str, str] = field(default_factory=dict)


@dataclass
class Job:
    """A tool and mode applied to one input (a file or a whole dataset)."""

    dataset: str
    scope: str
    input: str
    files: list[Path]
    tool: str
    mode: str
    threads: int
    commands: list[Command]

    @property
    def label(self) -> str:
        return f"{self.tool} {self.mode}"

    @property
    def input_mb(self) -> float:
        return sum(f.stat().st_size for f in self.files) / 1e6


class ReadSummaryBench:
    """Build and run read summary jobs."""

    def __init__(self, config: ReadSummaryConfig, out_dir: Path, segul_bin: Path):
        self.config = config
        self.out_dir = out_dir
        self.work_dir = out_dir / "work"
        self.segul_bin = segul_bin

    def jobs(
        self,
        datasets: list[Dataset],
        scopes: list[str],
        tools: list[str],
    ) -> list[list[Job]]:
        """Jobs grouped by input. Jobs in a group are interleaved per repeat."""
        groups = []
        for dataset in datasets:
            files = dataset.files()
            if not files:
                raise FileNotFoundError(
                    f"No FASTQ files for dataset {dataset.id} in {dataset.path}"
                )
            if "file" in scopes:
                for path in files:
                    name = path.name.removesuffix(".fq.gz")
                    groups.append(
                        self._group(dataset, "file", name, [path], [[path]], tools)
                    )
            if "dataset" in scopes:
                pairs = [list(p) for p in dataset.read_pairs()]
                groups.append(
                    self._group(dataset, "dataset", dataset.id, files, pairs, tools)
                )
        return groups

    def _group(
        self,
        dataset: Dataset,
        scope: str,
        name: str,
        files: list[Path],
        fastp_inputs: list[list[Path]],
        tools: list[str],
    ) -> list[Job]:
        base = self.work_dir / dataset.id / scope / name
        jobs = []
        threads = self.config.threads
        fastqc_threads = min(threads, FASTQC_THREADS_PER_FILE * len(files))
        if "segul" in tools:
            for mode in SEGUL_MODES:
                workdir = base / f"segul-{mode}"
                command = self._segul_command(files, mode, workdir, threads)
                jobs.append(
                    Job(
                        dataset.id,
                        scope,
                        name,
                        files,
                        "segul",
                        mode,
                        threads,
                        [command],
                    )
                )
            if "fastqc" in tools and fastqc_threads != threads:
                # Same thread budget as FastQC for the like-for-like comparison.
                workdir = base / f"segul-complete-{fastqc_threads}t"
                command = self._segul_command(
                    files, "complete", workdir, fastqc_threads
                )
                jobs.append(
                    Job(
                        dataset.id,
                        scope,
                        name,
                        files,
                        "segul",
                        "complete",
                        fastqc_threads,
                        [command],
                    )
                )
        if "fastp" in tools:
            for mode in FASTP_MODES:
                workdir = base / f"fastp-{mode}"
                commands = [
                    self._fastp_command(reads, mode, workdir) for reads in fastp_inputs
                ]
                jobs.append(
                    Job(
                        dataset.id, scope, name, files, "fastp", mode, threads, commands
                    )
                )
        if "fastqc" in tools:
            command = self._fastqc_command(
                files, base / "fastqc-complete", fastqc_threads
            )
            jobs.append(
                Job(
                    dataset.id,
                    scope,
                    name,
                    files,
                    "fastqc",
                    "complete",
                    fastqc_threads,
                    [command],
                )
            )
        return jobs

    def _segul_command(
        self, files: list[Path], mode: str, workdir: Path, threads: int
    ) -> Command:
        argv = [
            str(self.segul_bin),
            "read",
            "summary",
            "-i",
            *map(str, files),
            "--mode",
            mode,
            "--dup-size",
            str(self.config.dup_size),
            "-o",
            str(workdir / "out"),
            "--force",
        ]
        env = {"RAYON_NUM_THREADS": str(threads)}
        return Command(argv, workdir, env)

    def _fastp_command(self, reads: list[Path], mode: str, workdir: Path) -> Command:
        stem = reads[0].name.removesuffix(".fq.gz")
        argv = [str(self.config.fastp_bin), "-i", str(reads[0])]
        if len(reads) == 2:
            argv += ["-I", str(reads[1])]
        argv += [
            *FASTP_BASE_ARGS,
            *FASTP_MODE_ARGS[mode],
            "-w",
            str(self.config.threads),
            "-j",
            str(workdir / f"{stem}.json"),
            "-h",
            str(workdir / f"{stem}.html"),
        ]
        return Command(argv, workdir)

    def _fastqc_command(
        self, files: list[Path], workdir: Path, threads: int
    ) -> Command:
        argv = [
            str(self.config.fastqc_bin),
            "--threads",
            str(threads),
            "--limits",
            str(self.fastqc_limits()),
            "--dup_length",
            str(self.config.fastqc_dup_length),
            "--quiet",
            "--outdir",
            str(workdir),
            *map(str, files),
        ]
        return Command(argv, workdir)

    def fastqc_limits(self) -> Path:
        """FastQC limits file with the modules SEGUL lacks set to ignore.

        Built from FastQC's own default limits so warn/error thresholds stay
        unchanged.
        """
        path = self.work_dir / "fastqc-limits.txt"
        if path.exists():
            return path
        default = self.config.fastqc_bin.parent / "Configuration" / "limits.txt"
        lines = []
        for line in default.read_text().splitlines():
            parts = line.split()
            if len(parts) == 3 and parts[1] == "ignore":
                ignore = "1" if parts[0] in self.config.fastqc_ignore else "0"
                line = f"{parts[0]}\tignore\t{ignore}"
            lines.append(line)
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text("\n".join(lines) + "\n")
        return path

    def run_job(self, job: Job) -> tuple[Measurement, int | None]:
        """Run every command of a job and combine the measurements."""
        total = Measurement(0.0, 0.0, 0.0, 0.0, 0)
        reads: int | None = 0 if job.tool in ("fastp", "fastqc") else None
        for i, cmd in enumerate(job.commands):
            cmd.workdir.mkdir(parents=True, exist_ok=True)
            m = measure(cmd.argv, cmd.workdir, cmd.workdir / f"run-{i}", cmd.env)
            total.wall_s += m.wall_s
            total.user_s += m.user_s
            total.sys_s += m.sys_s
            total.max_rss_mb = max(total.max_rss_mb, m.max_rss_mb)
            total.exit_code = max(total.exit_code, m.exit_code)
            if job.tool == "fastp":
                reads += fastp_read_count(Path(cmd.argv[cmd.argv.index("-j") + 1]))
            elif job.tool == "fastqc" and m.exit_code == 0:
                reads += sum(fastqc_read_count(cmd.workdir, f) for f in job.files)
        return total, reads

    def run(self, groups: list[list[Job]], repeats: int, warmup: int, csv_path: Path):
        """Run all jobs, appending one CSV row per measured run."""
        n_runs = sum(len(g) for g in groups) * (repeats + warmup)
        done = 0
        with open(csv_path, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=FIELDS)
            writer.writeheader()
            for group in groups:
                for rep in range(-warmup, repeats):
                    for job in group:
                        done += 1
                        m, reads = self.run_job(job)
                        tag = "warmup" if rep < 0 else f"rep {rep + 1}/{repeats}"
                        print(
                            f"[{done}/{n_runs}] {job.dataset:10} {job.scope:7} "
                            f"{job.input[:40]:40} {job.label + f' {job.threads}t':19} {tag:10} "
                            f"{m.wall_s:8.2f} s {m.max_rss_mb:9.1f} MB"
                            + (f"  EXIT {m.exit_code}" if m.exit_code else ""),
                            flush=True,
                        )
                        if rep < 0:
                            continue
                        writer.writerow(
                            {
                                "dataset": job.dataset,
                                "scope": job.scope,
                                "input": job.input,
                                "n_files": len(job.files),
                                "input_mb": round(job.input_mb, 2),
                                "tool": job.tool,
                                "mode": job.mode,
                                "label": job.label,
                                "threads": job.threads,
                                "repeat": rep + 1,
                                "wall_s": round(m.wall_s, 4),
                                "user_s": m.user_s,
                                "sys_s": m.sys_s,
                                "max_rss_mb": round(m.max_rss_mb, 1),
                                "reads": "" if reads is None else reads,
                                "exit_code": m.exit_code,
                            }
                        )
                        f.flush()


def fastp_read_count(json_path: Path) -> int:
    """Total reads fastp saw before filtering."""
    report = json.loads(json_path.read_text())
    return report["summary"]["before_filtering"]["total_reads"]


def fastqc_read_count(workdir: Path, fastq: Path) -> int:
    """Total sequences from a FastQC report zip."""
    stem = fastq.name.removesuffix(".gz").removesuffix(".fq").removesuffix(".fastq")
    with zipfile.ZipFile(workdir / f"{stem}_fastqc.zip") as archive:
        data = archive.read(f"{stem}_fastqc/fastqc_data.txt").decode()
    for line in data.splitlines():
        if line.startswith("Total Sequences"):
            return int(line.split("\t")[1])
    raise ValueError(f"No read count in FastQC report for {fastq}")


def build_segul(repo: Path) -> Path:
    """Build the working tree in release mode and freeze a copy of the binary."""
    print(f"Building segul in {repo}", flush=True)
    subprocess.run(
        ["cargo", "build", "--release", "--bin", "segul"], cwd=repo, check=True
    )
    return freeze_segul(repo)


def freeze_segul(repo: Path) -> Path:
    """Copy the release binary to `.cache/bin` so later rebuilds don't affect a run."""
    built = repo / "target" / "release" / "segul"
    if not built.exists():
        raise FileNotFoundError(f"No release build at {built}")
    rev = git_rev(repo)
    dest = REPO_ROOT / ".cache" / "bin" / f"segul-{rev}"
    dest.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(built, dest)
    return dest


def git_rev(repo: Path) -> str:
    """Short commit hash, with `-dirty` when the working tree has changes."""
    rev = _output(["git", "rev-parse", "--short", "HEAD"], repo)
    dirty = _output(["git", "status", "--porcelain", "--untracked-files=no"], repo)
    return f"{rev}-dirty" if dirty else rev


def tool_version(argv: list[str]) -> str:
    """First line a tool prints for its version flag (stdout or stderr)."""
    proc = subprocess.run(argv, check=False, capture_output=True, text=True)
    return (proc.stdout or proc.stderr).strip().splitlines()[0]


def machine_info() -> dict[str, str | int]:
    """Host, CPU, memory, and OS of the benchmark machine."""
    cpu = ""
    for line in Path("/proc/cpuinfo").read_text().splitlines():
        if line.startswith("model name"):
            cpu = line.split(":", 1)[1].strip()
            break
    mem_kb = 0
    for line in Path("/proc/meminfo").read_text().splitlines():
        if line.startswith("MemTotal"):
            mem_kb = int(line.split()[1])
            break
    return {
        "host": platform.node(),
        "os": f"{platform.system()} {platform.release()}",
        "arch": platform.machine(),
        "cpu": cpu,
        "logical_cpus": os.cpu_count() or 0,
        "memory_gb": round(mem_kb / 1024**2, 1),
    }


def write_metadata(
    path: Path,
    config: ReadSummaryConfig,
    segul_bin: Path,
    groups: list[list[Job]],
    repeats: int,
    warmup: int,
):
    """Save machine, tool versions, settings, and example commands."""
    examples = {}
    for group in groups:
        for job in group:
            key = f"{job.scope}: {job.label} ({job.threads} threads)"
            examples.setdefault(key, " ".join(job.commands[0].argv))
    meta = {
        "date": datetime.now().isoformat(timespec="seconds"),
        "machine": machine_info(),
        "threads": config.threads,
        "repeats": repeats,
        "warmup": warmup,
        "segul": {
            "binary": str(segul_bin),
            "version": tool_version([str(segul_bin), "--version"]),
            "repo": str(config.segul_repo),
            "git_rev": git_rev(config.segul_repo),
            "dup_size": config.dup_size,
            "thread_setting": "RAYON_NUM_THREADS",
        },
        "fastp": {
            "binary": str(config.fastp_bin),
            "version": tool_version([str(config.fastp_bin), "--version"]),
            "thread_setting": "-w",
        },
        "fastqc": {
            "binary": str(config.fastqc_bin),
            "version": tool_version([str(config.fastqc_bin), "--version"]),
            "java": tool_version(["java", "-version"]),
            "thread_setting": "--threads (at most 4 per file)",
            "ignored_modules": config.fastqc_ignore,
            "dup_length": config.fastqc_dup_length,
        },
        "datasets": [
            {k: str(v) for k, v in asdict(d).items()}
            | {
                "n_files": len(d.files()),
                "size_mb": round(sum(f.stat().st_size for f in d.files()) / 1e6, 1),
            }
            for d in config.datasets
        ],
        "computes": config.computes,
        "example_commands": examples,
    }
    path.write_text(json.dumps(meta, indent=2) + "\n")


def _output(argv: list[str], cwd: Path) -> str:
    return subprocess.run(
        argv, cwd=cwd, capture_output=True, text=True, check=True
    ).stdout.strip()
