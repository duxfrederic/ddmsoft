import numpy as np
import pytest

from .fixtures import generate_model_data, read_legacy_matrix, write_legacy_matrix


def test_legacy_three_file_round_trip(tmp_path):
    source = generate_model_data("stretch")
    paths = write_legacy_matrix(tmp_path, "experiment", source)
    restored = read_legacy_matrix(tmp_path, "experiment")
    assert all(path.is_file() for path in paths)
    assert np.array_equal(restored.matrix, source.matrix)
    assert np.array_equal(restored.lag_times, source.lag_times)
    assert np.array_equal(restored.q_values, source.q_values)


def test_incomplete_legacy_matrix_is_reported(tmp_path):
    write_legacy_matrix(tmp_path, "experiment", generate_model_data("stretch"))
    (tmp_path / "ddm_matrices" / "experiment_QS.npy").unlink()
    with pytest.raises(FileNotFoundError, match="incomplete legacy matrix set"):
        read_legacy_matrix(tmp_path, "experiment")
