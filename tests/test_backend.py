import numpy as np
import pytest

from fitDDM import cumulant_exponential, fitOneDDMmatrix, merge
from generateDDM import FFTStack, logSpaced, partition_frame_counts
from utilities import RadialAverager, loadAnalyzedVideos, readParams
from fitDDM import mergeDDM


def test_log_spaced_supports_short_videos():
    for length in range(2, 10):
        lags = logSpaced(length, pointsPerDecade=10)
        assert len(lags) > 0
        assert lags.min() >= 1
        assert lags.max() < length


def test_log_spaced_rejects_invalid_input():
    with pytest.raises(ValueError):
        logSpaced(1)
    with pytest.raises(ValueError):
        logSpaced(10, pointsPerDecade=0)


def test_partition_counts_account_for_every_frame():
    assert partition_frame_counts(10, 3) == [4, 3, 3]
    with pytest.raises(ValueError):
        partition_frame_counts(3, 0)
    with pytest.raises(ValueError):
        partition_frame_counts(3, 4)


def test_radial_average_has_no_empty_bin_nans():
    averager = RadialAverager((8, 8), N=4)
    result = averager(np.ones((8, 8)))
    assert np.isfinite(result).all()


def test_load_video_from_generator_and_compute_matrix(tmp_path):
    frames = [np.full((8, 8), value, dtype=np.uint8) for value in range(4)]
    stack = FFTStack(freq=10, pixelsize=1e-6, maxCouples=0, ptPerDecade=10)
    stack.loadVideoFromGenerator(iter(frames), str(tmp_path / "sample.avi"), len(frames))
    stack.fftVideo()
    stack.stackToDDM()

    matrix_dir = tmp_path / "ddm_matrices"
    assert (matrix_dir / "sample_DDM_matrix.npy").exists()
    matrix = np.load(matrix_dir / "sample_DDM_matrix.npy")
    assert matrix.shape[0] == len(logSpaced(len(frames), 10))
    assert np.isfinite(matrix).all()


def test_read_params_ignores_comments_and_preserves_colons(tmp_path):
    config = tmp_path / "acquisition_parameters.txt"
    config.write_text(
        "# camera setup\n"
        "framerate: 120\n"
        "pixel note: camera: A\n"
        "\n",
        encoding="utf-8",
    )
    assert readParams(str(tmp_path)) == {
        "framerate": "120",
        "pixel note": "camera: A",
    }


def test_analyzed_matrix_display_names_map_to_full_paths(tmp_path):
    matrix_dir = tmp_path / "ddm_matrices"
    matrix_dir.mkdir()
    stem = matrix_dir / "sample"
    np.save(str(stem) + "_DDM_matrix.npy", np.ones((3, 2)))
    np.save(str(stem) + "_deltaTs.npy", np.array([1., 2., 3.]))
    np.save(str(stem) + "_QS.npy", np.array([1., 2.]))

    computed, displayed = loadAnalyzedVideos(str(tmp_path))

    full_path = str(stem) + "_DDM_matrix.npy"
    assert computed[full_path][0].shape == (3, 2)
    assert displayed["sample_DDM_matrix.npy"] == full_path


def test_average_matrices_groups_each_frame_rate(tmp_path):
    matrix_dir = tmp_path / "ddm_matrices"
    matrix_dir.mkdir()
    qs = np.array([1., 2.])
    datasets = {}
    for name, dts, value in (
        ("a", np.array([1., 2., 3.]), 1.),
        ("b", np.array([1., 2., 3.]), 3.),
        ("c", np.array([0.5, 1., 1.5]), 10.),
    ):
        path = str(matrix_dir / f"{name}_DDM_matrix.npy")
        datasets[path] = [np.full((3, 2), value), dts, qs]

    outputs = mergeDDM(datasets, mode="average")

    assert len(outputs) == 2
    matrix_outputs = [paths[0] for paths in outputs]
    assert any(np.allclose(np.load(path), 2.) for path in matrix_outputs)


def test_merge_does_not_mutate_inputs():
    dts_fast = np.array([1.0, 2.0, 3.0, 4.0])
    dts_slow = np.array([2.0, 4.0, 8.0])
    fast = np.arange(12.0).reshape(4, 3) + 1
    slow = np.arange(9.0).reshape(3, 3) + 2
    fast_before = fast.copy()
    slow_before = slow.copy()

    merged, dts = merge(fast, slow, dts_fast, dts_slow)

    assert np.array_equal(fast, fast_before)
    assert np.array_equal(slow, slow_before)
    assert np.all(np.diff(dts) > 0)
    assert merged.shape[1] == fast.shape[1]


def test_fit_returns_curve_matching_returned_parameters():
    dts = np.linspace(0.1, 2.0, 8)
    qs = np.array([1.0, 1.5, 2.0])
    params = np.array([2.0, 0.1, 0.7])
    QS, DTS = np.meshgrid(qs, dts)
    ddm = params[0] * (1 - cumulant_exponential(params, QS, DTS)) + params[1]

    amplitude, noise, cumulants, fitted, _ = fitOneDDMmatrix(
        (ddm, dts, qs),
        model="cumulant_1",
        ini=[0.7, "", ""],
        fixed=[False, False, False],
    )

    expected = cumulant_exponential(
        [amplitude[0], noise[0], cumulants[0][0]], qs[0], dts
    )
    assert np.allclose(fitted[:, 0], expected, rtol=1e-5, atol=1e-7)
