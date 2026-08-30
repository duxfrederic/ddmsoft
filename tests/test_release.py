from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

from ddmsoft.engine import VideoReadError, read_video_frames
from ddmsoft.io import MatrixLoadError, discover_matrix_sets, load_directory, load_matrices


def test_import_is_side_effect_free_and_resource_loads_outside_repository(tmp_path):
    script = """
import pathlib
import sys
import ddmsoft
from ddmsoft.resources import icon
assert ddmsoft.__version__ == '0.1.0'
assert icon().is_file()
assert not any(name.startswith('PySide6') for name in sys.modules)
assert 'tkinter' not in sys.modules
assert list(pathlib.Path.cwd().iterdir()) == []
"""

    completed = subprocess.run(
        [sys.executable, "-c", script],
        cwd=tmp_path,
        check=False,
        capture_output=True,
        text=True,
    )

    assert completed.returncode == 0, completed.stderr


def test_unicode_space_paths_and_read_only_inputs_remain_loadable(tmp_path):
    unicode_name = "\N{MICRO SIGN}"
    root = tmp_path / f"microscopy data {unicode_name}"
    root.mkdir()
    video = root / f"sample {unicode_name}.avi"
    metadata = root / f"sample {unicode_name}.txt"
    video.touch()
    metadata.write_text("framerate: 25\npixelsize: 1e-6\n", encoding="utf-8")
    matrix_directory = root / "ddm_matrices"
    matrix_directory.mkdir()
    arrays = (np.ones((2, 2)), np.array([0.1, 0.2]), np.array([1.0, 2.0]))
    paths = []
    for suffix, array in zip(("_DDM_matrix.npy", "_deltaTs.npy", "_QS.npy"), arrays):
        path = matrix_directory / f"sample {unicode_name}{suffix}"
        np.save(path, array)
        path.chmod(0o444)
        paths.append(path)
    video.chmod(0o444)
    metadata.chmod(0o444)

    assert load_directory(root)[video].frame_rate == 25
    loaded = load_matrices(root)[paths[0]]
    assert np.array_equal(loaded.matrix, arrays[0])


def test_legacy_directional_matrix_with_empty_bins_remains_loadable(tmp_path):
    root = tmp_path / "directional"
    matrix_directory = root / "ddm_matrices"
    matrix_directory.mkdir(parents=True)
    matrix = np.array([[1.0, np.nan], [2.0, 3.0]])
    np.save(matrix_directory / "sample_0.0__DDM_matrix.npy", matrix)
    np.save(matrix_directory / "sample_0.0__deltaTs.npy", np.array([0.1, 0.2]))
    np.save(matrix_directory / "sample_0.0__QS.npy", np.array([1.0, 2.0]))

    loaded = load_matrices(root)

    assert np.isnan(loaded[next(iter(loaded))].matrix[0, 1])


def test_corrupt_matrix_and_video_fail_descriptively(tmp_path):
    root = tmp_path / "corrupt"
    matrix_directory = root / "ddm_matrices"
    matrix_directory.mkdir(parents=True)
    np.save(matrix_directory / "bad_DDM_matrix.npy", np.ones(3))
    np.save(matrix_directory / "bad_deltaTs.npy", np.array([0.1, 0.2, 0.3]))
    np.save(matrix_directory / "bad_QS.npy", np.array([1.0]))

    matrix_set = discover_matrix_sets(root)[0]
    with pytest.raises(MatrixLoadError, match="could not load"):
        matrix_set.load()

    corrupt_video = tmp_path / "corrupt.avi"
    corrupt_video.write_bytes(b"not a video")
    with pytest.raises(VideoReadError, match="could not open|no decodable"):
        list(read_video_frames(corrupt_video))


def test_modern_source_declares_no_retired_gui_or_video_stack():
    root = Path(__file__).parents[1]
    source = "\n".join(path.read_text(encoding="utf-8") for path in (root / "src").rglob("*.py"))
    packaging = (root / "pyproject.toml").read_text(encoding="utf-8")
    retired = ("PySimpleGUI", "FreeSimpleGUI", "skvideo", "scikit-video", "tkinter", "TkAgg")

    assert all(name not in source for name in retired)
    assert all(name not in packaging for name in retired)
