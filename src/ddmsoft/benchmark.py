"""Run deterministic generated benchmarks against the modern DDM backends."""

from __future__ import annotations

import argparse
import gc
import json
import platform
import tempfile
import time
import tracemalloc
from collections.abc import Callable, Sequence
from numbers import Integral
from pathlib import Path
from typing import Any

import numpy as np

from .engine import compute_ddm
from .models import DDMData
from .time_dependent import compute_time_dependent_ddm

DEFAULT_FRAME_COUNT = 24
DEFAULT_FRAME_SIZE = 32
DEFAULT_SECTORS = 2
DEFAULT_PARTITIONS = 3
DEFAULT_MAX_COUPLES = 24
POINTS_PER_DECADE = 4
SEED = 15


def _positive_integer(value: int, name: str) -> int:
    if isinstance(value, bool) or not isinstance(value, Integral) or value <= 0:
        raise ValueError(f"{name} must be a positive integer")
    return int(value)


def _generated_frames(frame_count: int, frame_size: int) -> np.ndarray:
    rng = np.random.default_rng(SEED)
    texture = rng.normal(0.0, 1.0, (frame_size, frame_size))
    frames = [
        np.roll(texture, shift=(index // 2, index), axis=(0, 1))
        + 0.1 * rng.normal(0.0, 1.0, texture.shape)
        for index in range(frame_count)
    ]
    return np.asarray(frames, dtype=float)


def _measure(operation: Callable[[], Any]) -> tuple[Any, float, int]:
    gc.collect()
    tracemalloc.start()
    started = time.perf_counter()
    try:
        result = operation()
        runtime = time.perf_counter() - started
        _, peak = tracemalloc.get_traced_memory()
    finally:
        tracemalloc.stop()
    return result, runtime, peak


def _record(
    name: str,
    frame_count: int,
    frame_size: int,
    lag_count: int,
    sectors: int,
    partitions: int,
    runtime: float,
    peak: int,
) -> dict[str, object]:
    return {
        "frame_count": frame_count,
        "frame_dimensions": [frame_size, frame_size],
        "lag_count": lag_count,
        "name": name,
        "partitions": partitions,
        "peak_traced_memory_bytes": peak,
        "runtime_seconds": round(runtime, 9),
        "sectors": sectors,
    }


def run_benchmarks(
    *,
    frame_count: int = DEFAULT_FRAME_COUNT,
    frame_size: int = DEFAULT_FRAME_SIZE,
    sectors: int = DEFAULT_SECTORS,
    partitions: int = DEFAULT_PARTITIONS,
    max_couples: int = DEFAULT_MAX_COUPLES,
) -> dict[str, object]:
    """Run isotropic, directional, and time-dependent generated workloads."""
    frame_count = _positive_integer(frame_count, "frame_count")
    frame_size = _positive_integer(frame_size, "frame_size")
    sectors = _positive_integer(sectors, "sectors")
    partitions = _positive_integer(partitions, "partitions")
    if frame_count < 2:
        raise ValueError("frame_count must be at least two")
    if frame_size < 2:
        raise ValueError("frame_size must be at least two")
    if sectors < 2:
        raise ValueError("sectors must be at least two for the directional benchmark")
    if frame_count < 2 * partitions:
        raise ValueError("frame_count must provide at least two frames per partition")
    if isinstance(max_couples, bool) or not isinstance(max_couples, Integral) or max_couples < 0:
        raise ValueError("max_couples must be a non-negative integer")
    max_couples = int(max_couples)
    frames = _generated_frames(frame_count, frame_size)
    common = {
        "max_couples": max_couples,
        "points_per_decade": POINTS_PER_DECADE,
    }

    isotropic, runtime, peak = _measure(lambda: compute_ddm(frames, 20.0, 1.0e-6, **common))
    if not isinstance(isotropic, DDMData):
        raise TypeError("isotropic benchmark returned directional data")
    records = [
        _record(
            "isotropic",
            frame_count,
            frame_size,
            isotropic.matrix.shape[0],
            1,
            1,
            runtime,
            peak,
        )
    ]

    directional, runtime, peak = _measure(
        lambda: compute_ddm(frames, 20.0, 1.0e-6, sectors=sectors, **common)
    )
    if not isinstance(directional, tuple):
        raise TypeError("directional benchmark returned isotropic data")
    records.append(
        _record(
            "directional",
            frame_count,
            frame_size,
            directional[0].matrix.shape[0],
            len(directional),
            1,
            runtime,
            peak,
        )
    )

    time_dependent, runtime, peak = _measure(
        lambda: compute_time_dependent_ddm(
            frames,
            20.0,
            1.0e-6,
            partitions,
            **common,
        )
    )
    partition_lag_counts = [
        result.data.matrix.shape[0] for result in time_dependent if isinstance(result.data, DDMData)
    ]
    if len(partition_lag_counts) != partitions:
        raise RuntimeError("time-dependent benchmark returned directional data")
    time_record = _record(
        "time_dependent",
        frame_count,
        frame_size,
        sum(partition_lag_counts),
        1,
        len(time_dependent),
        runtime,
        peak,
    )
    time_record["lag_counts_per_partition"] = partition_lag_counts
    records.append(time_record)

    return {
        "configuration": {
            "max_couples": max_couples,
            "points_per_decade": POINTS_PER_DECADE,
            "seed": SEED,
        },
        "environment": {
            "memory_measurement": "tracemalloc peak traced bytes",
            "platform": platform.platform(),
            "python": platform.python_version(),
        },
        "schema_version": 1,
        "workloads": records,
    }


def write_results(path: str | Path, results: dict[str, object]) -> Path:
    """Atomically write benchmark results as consistently formatted JSON."""
    destination = Path(path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary_path: Path | None = None
    try:
        with tempfile.NamedTemporaryFile(
            mode="w",
            encoding="utf-8",
            newline="\n",
            prefix=f".{destination.name}.",
            suffix=".tmp",
            dir=destination.parent,
            delete=False,
        ) as temporary:
            temporary_path = Path(temporary.name)
            json.dump(results, temporary, indent=2, sort_keys=True)
            temporary.write("\n")
        temporary_path.replace(destination)
    finally:
        if temporary_path is not None:
            temporary_path.unlink(missing_ok=True)
    return destination


def _argument_integer(value: str) -> int:
    try:
        return _positive_integer(int(value), "value")
    except ValueError as error:
        raise argparse.ArgumentTypeError("must be a positive integer") from error


def _non_negative_integer(value: str) -> int:
    try:
        result = int(value)
    except ValueError as error:
        raise argparse.ArgumentTypeError("must be a non-negative integer") from error
    if result < 0:
        raise argparse.ArgumentTypeError("must be a non-negative integer")
    return result


def main(argv: Sequence[str] | None = None) -> int:
    """Run generated benchmarks and write their JSON measurements."""
    parser = argparse.ArgumentParser(
        prog="ddmsoft-benchmark",
        description="Benchmark deterministic generated DDM workloads.",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        default=Path("ddmsoft-benchmark.json"),
        help="JSON output path (default: %(default)s)",
    )
    parser.add_argument(
        "--frame-count",
        "--frames",
        type=_argument_integer,
        default=DEFAULT_FRAME_COUNT,
        help="generated source frame count (default: %(default)s)",
    )
    parser.add_argument(
        "--frame-size",
        "--size",
        type=_argument_integer,
        default=DEFAULT_FRAME_SIZE,
        help="square frame width and height (default: %(default)s)",
    )
    parser.add_argument(
        "--sectors",
        type=_argument_integer,
        default=DEFAULT_SECTORS,
        help="directional sector count (default: %(default)s)",
    )
    parser.add_argument(
        "--partitions",
        type=_argument_integer,
        default=DEFAULT_PARTITIONS,
        help="time-dependent partition count (default: %(default)s)",
    )
    parser.add_argument(
        "--max-couples",
        type=_non_negative_integer,
        default=DEFAULT_MAX_COUPLES,
        help="maximum frame couples, or zero for all (default: %(default)s)",
    )
    arguments = parser.parse_args(argv)
    try:
        results = run_benchmarks(
            frame_count=arguments.frame_count,
            frame_size=arguments.frame_size,
            sectors=arguments.sectors,
            partitions=arguments.partitions,
            max_couples=arguments.max_couples,
        )
        output = write_results(arguments.output, results)
    except (OSError, RuntimeError, ValueError) as error:
        parser.error(str(error))
    print(output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
