from __future__ import annotations

import numpy as np
import pytest

from ddmsoft.io import (
    AmbiguousMetadataError,
    IncompleteMatrixError,
    InvalidMetadataError,
    MetadataFileError,
    MissingMetadataError,
    discover_matrix_sets,
    load_directory,
    load_matrices,
    parse_metadata_file,
    save_autocorrelation_csv,
    save_fit_text,
    save_matrix_csv,
)
from ddmsoft.models import DDMData


def _write(path, text: str) -> None:
    path.write_text(text, encoding="utf-8")


def test_metadata_ignores_comments_and_splits_only_first_colon(tmp_path):
    metadata = tmp_path / "acquisition_parameters.txt"
    _write(metadata, "\n# ignored\nframerate: 25\npixelsize: 1.5\nnote: room: A\n")
    assert parse_metadata_file(metadata) == {
        "framerate": "25",
        "pixelsize": "1.5",
        "note": "room: A",
    }
    (tmp_path / "sample.avi").touch()
    loaded = load_directory(tmp_path)
    assert loaded[tmp_path / "sample.avi"].frame_rate == 25
    assert loaded[tmp_path / "sample.avi"].pixel_size == 1.5


def test_metadata_supports_one_file_per_video_and_rejects_ambiguous_files(tmp_path):
    for name in ("a.avi", "b.avi"):
        (tmp_path / name).touch()
    _write(tmp_path / "a.txt", "framerate: 10\npixelsize: 1\n")
    _write(tmp_path / "b.txt", "framerate: 20\npixelsize: 2\n")
    loaded = load_directory(tmp_path)
    assert loaded[tmp_path / "a.avi"].frame_rate == 10
    assert loaded[tmp_path / "b.avi"].pixel_size == 2
    _write(tmp_path / "extra.txt", "framerate: 1\npixelsize: 1\n")
    with pytest.raises(AmbiguousMetadataError):
        load_directory(tmp_path)


@pytest.mark.parametrize(
    "text, error",
    [
        ("framerate 10\npixelsize: 1\n", MetadataFileError),
        ("framerate: nope\npixelsize: 1\n", InvalidMetadataError),
        ("framerate: 10\n", MissingMetadataError),
    ],
)
def test_metadata_values_are_validated_separately(tmp_path, text, error):
    (tmp_path / "sample.avi").touch()
    _write(tmp_path / "params.txt", text)
    with pytest.raises(error):
        load_directory(tmp_path)


def test_matrix_discovery_preserves_legacy_names_and_reports_incomplete_sets(tmp_path):
    matrix_dir = tmp_path / "ddm_matrices"
    matrix_dir.mkdir()
    arrays = (np.ones((2, 3)), np.array([0.1, 0.2]), np.array([1.0, 2.0, 3.0]))
    for suffix, array in zip(("_DDM_matrix.npy", "_deltaTs.npy", "_QS.npy"), arrays):
        np.save(matrix_dir / f"run_0.0_{suffix}", array)
    data = load_matrices(tmp_path)[matrix_dir / "run_0.0__DDM_matrix.npy"]
    assert np.array_equal(data.matrix, arrays[0])
    (matrix_dir / "broken_DDM_matrix.npy").touch()
    with pytest.raises(IncompleteMatrixError):
        discover_matrix_sets(tmp_path)
    assert len(discover_matrix_sets(tmp_path, strict=False)) == 1


def test_exports_use_exact_suffixes_and_utf8_text(tmp_path):
    data = DDMData(np.array([[1.0, 2.0], [3.0, 4.0]]), np.array([0.1, 0.2]), np.array([1.0, 2.0]))
    paths = save_matrix_csv(tmp_path / "result", data)
    assert [path.name for path in paths] == [
        "result_DDM_matrix.csv",
        "result_deltaTs.csv",
        "result_QS.csv",
    ]
    assert save_matrix_csv(tmp_path / "result_DDM_matrix.csv", data)[0] == paths[0]
    correlation = save_autocorrelation_csv(
        tmp_path / "corr", data, np.array([2.0, 3.0]), np.array([0.1, 0.2]), [2.0]
    )
    assert [path.name for path in correlation] == [
        "corr_autocorrelationmatrix.csv",
        "corr_qs.csv",
        "corr_dts.csv",
    ]
    fit_path = save_fit_text(tmp_path / "fit", [1.0], [2.0], [0.1], [np.array([1.0e-12])], ["D"])
    assert fit_path.name == "fit.txt"
    assert fit_path.read_text(encoding="utf-8").splitlines()[0] == "q [m^-1]\tA\tB\tD"
