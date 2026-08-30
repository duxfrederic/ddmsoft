import numpy as np

from .fixtures import (
    FIT_MODEL_IDS,
    changing_intensity_frames,
    constant_frames,
    direct_ddm_reference,
    generate_model_data,
    orthogonal_gratings,
    random_frames,
    translated_sinusoidal_grating,
)


def test_frame_generators_are_deterministic():
    assert np.array_equal(random_frames(seed=7), random_frames(seed=7))
    assert np.all(direct_ddm_reference(constant_frames()) == 0)
    assert changing_intensity_frames().shape == (5, 8, 8)
    assert translated_sinusoidal_grating().shape == (6, 16, 16)
    assert orthogonal_gratings().shape == (6, 16, 16)


def test_direct_reference_matches_its_definition():
    frames = random_frames(count=4, shape=(4, 5), seed=3)
    actual = direct_ddm_reference(frames, [1, 3])
    transformed = np.fft.fft2(frames, axes=(1, 2))
    expected = np.asarray(
        [
            np.mean(np.abs(transformed[lag:] - transformed[:-lag]) ** 2, axis=0)
            for lag in (1, 3)
        ]
    )
    assert np.allclose(actual, expected)


def test_every_legacy_fit_model_has_a_generated_matrix():
    generated = {model_id: generate_model_data(model_id, noise=1.0e-5) for model_id in FIT_MODEL_IDS}
    assert set(generated) == set(FIT_MODEL_IDS)
    assert all(data.matrix.shape == (7, 4) for data in generated.values())
    assert np.array_equal(
        generate_model_data("cumulant_1", noise=1.0e-5).matrix,
        generate_model_data("cumulant_1", noise=1.0e-5).matrix,
    )
