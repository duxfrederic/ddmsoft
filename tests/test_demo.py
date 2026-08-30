from __future__ import annotations

import numpy as np
import pytest

from ddmsoft.demo import FRAME_COUNT, FRAME_SIZE, generate_demo, generate_frames, main
from ddmsoft.engine import read_video_frames
from ddmsoft.io import discover_matrix_sets, load_directory


def test_demo_is_deterministic_valid_and_preserves_unrelated_files(tmp_path):
    output = tmp_path / "demo output"
    output.mkdir()
    unrelated = output / "keep-me.dat"
    unrelated.write_bytes(b"unrelated")

    with pytest.raises(FileExistsError, match="not empty"):
        generate_demo(output)

    assert generate_demo(output, overwrite=True) == output
    assert unrelated.read_bytes() == b"unrelated"
    assert np.array_equal(generate_frames(), generate_frames())

    metadata = load_directory(output)
    video = output / "demo.avi"
    assert metadata[video].frame_rate == 20.0
    assert metadata[video].pixel_size == 1.0e-6
    decoded = list(read_video_frames(video))
    assert len(decoded) == FRAME_COUNT
    assert all(frame.shape == (FRAME_SIZE, FRAME_SIZE) for frame in decoded)

    matrix_sets = discover_matrix_sets(output)
    assert len(matrix_sets) == 1
    data = matrix_sets[0].load()
    assert data.matrix.shape == (data.lag_times.size, data.q_values.size)
    assert np.all(np.isfinite(data.matrix))


def test_demo_cli_requires_overwrite_for_nonempty_output(tmp_path):
    output = tmp_path / "demo"
    output.mkdir()
    (output / "existing.txt").touch()

    with pytest.raises(SystemExit) as error:
        main(["--output", str(output)])
    assert error.value.code == 2
