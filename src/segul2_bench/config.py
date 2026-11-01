"""Load benchmark configs from `bench/*.toml`."""

import tomllib
from dataclasses import dataclass, field
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]


def resolve_path(value: str) -> Path:
    """Expand `~` and make relative paths relative to the repo root."""
    path = Path(value).expanduser()
    return path if path.is_absolute() else REPO_ROOT / path


@dataclass
class Dataset:
    """A folder of gzipped FASTQ files."""

    id: str
    description: str
    path: Path
    glob: str
    pairs: tuple[str, str]

    def files(self) -> list[Path]:
        """FASTQ files in the dataset, sorted by path."""
        return sorted(self.path.glob(self.glob))

    def read_pairs(self) -> list[tuple[Path, Path]]:
        """R1/R2 file pairs, matched by the configured suffixes."""
        r1_suffix, r2_suffix = self.pairs
        pairs = []
        for r1 in self.files():
            if r1.name.endswith(r1_suffix):
                r2 = r1.with_name(r1.name.removesuffix(r1_suffix) + r2_suffix)
                if not r2.exists():
                    raise FileNotFoundError(f"Missing read 2 for {r1}")
                pairs.append((r1, r2))
        return pairs


@dataclass
class ReadSummaryConfig:
    """Settings for the read summary benchmark."""

    threads: int
    repeats: int
    warmup: int
    segul_repo: Path
    dup_size: int
    fastp_bin: Path
    fastqc_bin: Path
    fastqc_ignore: list[str]
    fastqc_dup_length: int
    datasets: list[Dataset]
    computes: dict[str, str] = field(default_factory=dict)


def load_read_summary(path: Path) -> ReadSummaryConfig:
    """Read `bench/read-summary.toml`."""
    with open(path, "rb") as f:
        raw = tomllib.load(f)
    defaults = raw["defaults"]
    tools = raw["tools"]
    datasets = [
        Dataset(
            id=d["id"],
            description=d.get("description", ""),
            path=resolve_path(d["path"]),
            glob=d["glob"],
            pairs=tuple(d["pairs"]),
        )
        for d in raw["datasets"]
    ]
    return ReadSummaryConfig(
        threads=defaults["threads"],
        repeats=defaults["repeats"],
        warmup=defaults.get("warmup", 1),
        segul_repo=resolve_path(tools["segul"]["segul_repo"]),
        dup_size=tools["segul"].get("dup_size", 100_000),
        fastp_bin=resolve_path(tools["fastp"]["binary"]),
        fastqc_bin=resolve_path(tools["fastqc"]["binary"]),
        fastqc_ignore=tools["fastqc"].get("ignore", []),
        fastqc_dup_length=tools["fastqc"].get("dup_length", 0),
        datasets=datasets,
        computes=raw.get("computes", {}),
    )
