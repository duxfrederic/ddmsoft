"""Pure operations for combining computed DDM datasets."""

from __future__ import annotations

from collections.abc import Iterable, Sequence
from dataclasses import dataclass

import numpy as np

from .models import DDMData


class CombinationError(ValueError):
    """Base class for invalid matrix-combination inputs."""


@dataclass(frozen=True)
class NamedDDMData:
    """A deterministic display name paired with a computed dataset."""

    name: str
    data: DDMData


def _as_sequence(datasets: Iterable[DDMData]) -> tuple[DDMData, ...]:
    values = tuple(datasets)
    if not values:
        raise CombinationError("at least one DDM dataset is required")
    if not all(isinstance(value, DDMData) for value in values):
        raise TypeError("all datasets must be DDMData instances")
    return values


def _validate_monotonic(data: DDMData) -> None:
    if data.lag_times.size < 1 or np.any(np.diff(data.lag_times) <= 0):
        raise CombinationError(
            f"lag times must be finite and strictly increasing for {data.source or 'dataset'}"
        )


def _validate_compatible(datasets: Sequence[DDMData]) -> None:
    first = datasets[0]
    _validate_monotonic(first)
    for index, data in enumerate(datasets[1:], start=1):
        _validate_monotonic(data)
        if data.matrix.ndim != 2 or data.matrix.shape[1] != first.matrix.shape[1]:
            raise CombinationError(
                f"dataset {index} has incompatible matrix width: "
                f"expected {first.matrix.shape[1]}, got {data.matrix.shape[1]}"
            )
        if not np.allclose(data.q_values, first.q_values, rtol=1e-9, atol=0.0):
            raise CombinationError(f"dataset {index} has an incompatible q grid")


def _merge_pair(fast: DDMData, slow: DDMData) -> DDMData:
    """Merge two ordered datasets using the legacy scale-and-transition rule."""
    _validate_compatible((fast, slow))
    if fast.lag_times[0] > slow.lag_times[0]:
        fast, slow = slow, fast
    if fast.lag_times.shape == slow.lag_times.shape and np.allclose(
        fast.lag_times, slow.lag_times, rtol=1e-9, atol=0.0
    ):
        raise CombinationError("cannot merge datasets with identical lag-time grids")

    fast_times = fast.lag_times
    slow_times = slow.lag_times
    overlap_mask = (fast_times >= slow_times[0]) & (fast_times <= slow_times[-1])
    if not np.any(overlap_mask):
        if fast_times[-1] < slow_times[0]:
            return DDMData(
                np.vstack((fast.matrix, slow.matrix)),
                np.concatenate((fast_times, slow_times)),
                fast.q_values,
            )
        if slow_times[-1] < fast_times[0]:
            return DDMData(
                np.vstack((slow.matrix, fast.matrix)),
                np.concatenate((slow_times, fast_times)),
                fast.q_values,
            )
        raise CombinationError("lag-time grids overlap without an ordered transition")

    first_overlap = int(np.flatnonzero(overlap_mask)[0])
    overlap_times = fast_times[overlap_mask]
    slow_first = slow.matrix[0]
    fast_at_start = fast.matrix[first_overlap]
    zero_slow = np.isclose(slow_first, 0.0)
    zero_fast = np.isclose(fast_at_start, 0.0)
    if np.any(zero_slow & ~zero_fast):
        raise CombinationError(
            "cannot scale the slow matrix because its first row contains zero values"
        )
    scale = np.ones_like(slow_first, dtype=float)
    np.divide(fast_at_start, slow_first, out=scale, where=~zero_slow)
    scaled_slow = slow.matrix * scale[None, :]
    interpolated = np.column_stack(
        [
            np.interp(overlap_times, slow_times, scaled_slow[:, column])
            for column in range(slow.matrix.shape[1])
        ]
    )
    if overlap_times.size == 1:
        transition = interpolated
    else:
        blend = np.linspace(0.0, 1.0, overlap_times.size)[:, None]
        transition = (1.0 - blend) * fast.matrix[overlap_mask] + blend * interpolated
    slow_after = slow_times > overlap_times[-1]
    merged_times = np.concatenate(
        (fast_times[:first_overlap], overlap_times, slow_times[slow_after])
    )
    merged_matrix = np.vstack(
        (fast.matrix[:first_overlap], transition, scaled_slow[slow_after])
    )
    return DDMData(merged_matrix, merged_times, fast.q_values)


def merge_ddm(datasets: Iterable[DDMData]) -> DDMData:
    """Merge datasets without mutating any input arrays."""
    values = _as_sequence(datasets)
    _validate_compatible(values)
    ordered = sorted(values, key=lambda value: value.lag_times[0])
    result = ordered[0]
    for data in ordered[1:]:
        result = _merge_pair(result, data)
    return DDMData(result.matrix, result.lag_times, result.q_values)


def average_ddm(datasets: Iterable[DDMData]) -> DDMData:
    """Average compatible datasets without modifying their arrays."""
    values = _as_sequence(datasets)
    _validate_compatible(values)
    first = values[0]
    if any(
        data.matrix.shape != first.matrix.shape
        or not np.allclose(data.lag_times, first.lag_times, rtol=1e-9, atol=0.0)
        for data in values[1:]
    ):
        raise CombinationError(
            "datasets must have compatible matrix shapes and lag-time grids to average"
        )
    matrix = np.mean(np.stack([data.matrix for data in values], axis=0), axis=0)
    return DDMData(matrix, first.lag_times, first.q_values)


def average_groups(datasets: Iterable[DDMData]) -> tuple[NamedDDMData, ...]:
    """Group by compatible lag arrays and return one named result per group."""
    values = _as_sequence(datasets)
    _validate_compatible(values)
    groups: list[list[DDMData]] = []
    for data in values:
        for group in groups:
            if (
                data.matrix.shape == group[0].matrix.shape
                and data.lag_times.shape == group[0].lag_times.shape
            ):
                same_lags = np.allclose(
                    data.lag_times, group[0].lag_times, rtol=1e-9, atol=0.0
                )
            else:
                same_lags = False
            if same_lags:
                group.append(data)
                break
        else:
            groups.append([data])
    return tuple(
        NamedDDMData(f"average_{index}", average_ddm(group))
        for index, group in enumerate(groups)
    )


# Compact aliases for callers migrating from the legacy terminology.
merge_data = merge_ddm
average_data = average_ddm
group_averages = average_groups


__all__ = [
    "CombinationError",
    "NamedDDMData",
    "average_data",
    "average_ddm",
    "average_groups",
    "group_averages",
    "merge_data",
    "merge_ddm",
]
