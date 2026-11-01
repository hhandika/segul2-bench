"""Run a command under `/usr/bin/time -v` and collect time and memory."""

import os
import subprocess
import time
from dataclasses import dataclass
from pathlib import Path

TIME_BIN = Path("/usr/bin/time")


@dataclass
class Measurement:
    """Resources used by one command."""

    wall_s: float
    user_s: float
    sys_s: float
    max_rss_mb: float
    exit_code: int


def measure(
    cmd: list[str],
    cwd: Path,
    log_prefix: Path,
    env: dict[str, str] | None = None,
) -> Measurement:
    """Run `cmd` in `cwd`, writing its stdout/stderr next to `log_prefix`.

    Wall time comes from `time.perf_counter` (finer than GNU time's 10 ms
    resolution). CPU time and peak RSS come from GNU time.
    """
    time_out = log_prefix.with_suffix(".time")
    full_env = {**os.environ, **(env or {})}
    with (
        open(log_prefix.with_suffix(".stdout"), "w") as out,
        open(log_prefix.with_suffix(".stderr"), "w") as err,
    ):
        start = time.perf_counter()
        proc = subprocess.run(
            [str(TIME_BIN), "-v", "-o", str(time_out), *cmd],
            check=False,
            cwd=cwd,
            env=full_env,
            stdout=out,
            stderr=err,
            stdin=subprocess.DEVNULL,
        )
        wall = time.perf_counter() - start
    stats = parse_time_output(time_out.read_text())
    return Measurement(
        wall_s=wall,
        user_s=stats["user_s"],
        sys_s=stats["sys_s"],
        max_rss_mb=stats["max_rss_kb"] / 1024,
        exit_code=proc.returncode,
    )


def parse_time_output(text: str) -> dict[str, float]:
    """Parse the fields we need from `/usr/bin/time -v` output."""
    keys = {
        "User time (seconds)": "user_s",
        "System time (seconds)": "sys_s",
        "Maximum resident set size (kbytes)": "max_rss_kb",
    }
    stats = {value: 0.0 for value in keys.values()}
    for line in text.splitlines():
        name, _, value = line.strip().rpartition(": ")
        if name in keys:
            stats[keys[name]] = float(value)
    return stats
