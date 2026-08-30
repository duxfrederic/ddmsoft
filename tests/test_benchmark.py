from __future__ import annotations

import json

import pytest

from ddmsoft.benchmark import main, run_benchmarks, write_results


def test_generated_benchmarks_record_all_required_measurements(tmp_path):
    results = run_benchmarks(
        frame_count=8,
        frame_size=8,
        sectors=2,
        partitions=2,
        max_couples=4,
    )

    assert results["schema_version"] == 1
    workloads = results["workloads"]
    assert [record["name"] for record in workloads] == [
        "isotropic",
        "directional",
        "time_dependent",
    ]
    assert [record["sectors"] for record in workloads] == [1, 2, 1]
    assert [record["partitions"] for record in workloads] == [1, 1, 2]
    for record in workloads:
        assert record["frame_count"] == 8
        assert record["frame_dimensions"] == [8, 8]
        assert record["lag_count"] > 0
        assert record["runtime_seconds"] >= 0
        assert record["peak_traced_memory_bytes"] > 0

    output = write_results(tmp_path / "results" / "benchmark.json", results)
    assert json.loads(output.read_text(encoding="utf-8")) == results
    assert output.read_bytes().endswith(b"\n")


def test_benchmark_cli_accepts_output_and_size_controls(monkeypatch, tmp_path):
    captured = {}

    def fake_run(**kwargs):
        captured.update(kwargs)
        return {"schema_version": 1, "workloads": []}

    monkeypatch.setattr("ddmsoft.benchmark.run_benchmarks", fake_run)
    output = tmp_path / "benchmark.json"
    assert (
        main(
            [
                "--output",
                str(output),
                "--frames",
                "10",
                "--size",
                "12",
                "--sectors",
                "3",
                "--partitions",
                "2",
                "--max-couples",
                "0",
            ]
        )
        == 0
    )
    assert captured == {
        "frame_count": 10,
        "frame_size": 12,
        "sectors": 3,
        "partitions": 2,
        "max_couples": 0,
    }
    assert json.loads(output.read_text(encoding="utf-8"))["schema_version"] == 1


def test_benchmark_rejects_partitions_with_fewer_than_two_frames():
    with pytest.raises(ValueError, match="two frames per partition"):
        run_benchmarks(frame_count=5, frame_size=8, partitions=3)

    with pytest.raises(ValueError, match="at least two"):
        run_benchmarks(frame_count=8, frame_size=8, sectors=1, partitions=2)
