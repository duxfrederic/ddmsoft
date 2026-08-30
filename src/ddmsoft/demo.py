"""Generate a small, self-contained DDMSoft demonstration dataset."""

from __future__ import annotations

import argparse
from collections.abc import Sequence
from pathlib import Path
from tempfile import TemporaryDirectory

import cv2
import numpy as np

from .engine import compute_video_ddm
from .io import DDM_MATRICES_DIRECTORY, save_matrix_set
from .models import DDMData

FRAME_COUNT = 12
FRAME_RATE = 20.0
FRAME_SIZE = 32
PIXEL_SIZE = 1.0e-6
VIDEO_NAME = "demo.avi"
METADATA_NAME = "demo.txt"
MATRIX_STEM = "demo"


def generate_frames() -> np.ndarray:
    """Return the deterministic grayscale frames used by the demo video."""
    y, x = np.indices((FRAME_SIZE, FRAME_SIZE), dtype=float)
    frames = []
    for index in range(FRAME_COUNT):
        horizontal = np.sin(2.0 * np.pi * (x - index) / FRAME_SIZE)
        vertical = np.cos(4.0 * np.pi * (y - 2 * index) / FRAME_SIZE)
        spot_x = (3 * index + 7) % FRAME_SIZE
        spot_y = (2 * index + 11) % FRAME_SIZE
        distance = (x - spot_x) ** 2 + (y - spot_y) ** 2
        image = 112.0 + 42.0 * horizontal + 30.0 * vertical + 55.0 * np.exp(-distance / 18.0)
        frames.append(np.clip(np.rint(image), 0, 255).astype(np.uint8))
    return np.asarray(frames)


def _video_is_valid(path: Path, expected_frames: int) -> bool:
    capture = cv2.VideoCapture(str(path))
    decoded = 0
    try:
        if not capture.isOpened():
            return False
        while True:
            ok, frame = capture.read()
            if not ok:
                break
            if frame.shape[:2] != (FRAME_SIZE, FRAME_SIZE):
                return False
            decoded += 1
    finally:
        capture.release()
    return decoded == expected_frames


def _write_video(path: Path, frames: np.ndarray) -> None:
    for codec in ("MJPG", "XVID", "mp4v"):
        writer = cv2.VideoWriter(
            str(path),
            cv2.VideoWriter_fourcc(*codec),
            FRAME_RATE,
            (FRAME_SIZE, FRAME_SIZE),
            True,
        )
        try:
            if not writer.isOpened():
                continue
            for frame in frames:
                writer.write(cv2.cvtColor(frame, cv2.COLOR_GRAY2BGR))
        except cv2.error:
            continue
        finally:
            writer.release()
        if _video_is_valid(path, len(frames)):
            return
        path.unlink(missing_ok=True)
    raise RuntimeError("OpenCV could not create a decodable AVI with an available codec")


def generate_demo(output_directory: str | Path, *, overwrite: bool = False) -> Path:
    """Generate the AVI, metadata, and legacy matrix set below one directory."""
    output = Path(output_directory)
    if output.exists():
        if not output.is_dir():
            raise FileExistsError(f"output path is not a directory: {output}")
        if any(output.iterdir()) and not overwrite:
            raise FileExistsError(f"output directory is not empty: {output}")

    output.parent.mkdir(parents=True, exist_ok=True)
    frames = generate_frames()
    with TemporaryDirectory(prefix=".ddmsoft-demo-", dir=output.parent) as temporary:
        staging = Path(temporary)
        video = staging / VIDEO_NAME
        _write_video(video, frames)
        (staging / METADATA_NAME).write_text(
            "# Generated DDMSoft demonstration acquisition\n"
            f"framerate: {FRAME_RATE:g}\n"
            f"pixelsize: {PIXEL_SIZE:g}\n"
            "temperature: 25\n",
            encoding="utf-8",
            newline="\n",
        )
        matrix_directory = staging / DDM_MATRICES_DIRECTORY
        matrix_directory.mkdir()
        data = compute_video_ddm(
            video,
            FRAME_RATE,
            PIXEL_SIZE,
            max_couples=0,
            points_per_decade=4,
        )
        if not isinstance(data, DDMData):
            raise TypeError("isotropic demo computation returned directional data")
        matrix_paths = save_matrix_set(matrix_directory / MATRIX_STEM, data)

        output.mkdir(parents=True, exist_ok=True)
        destination_matrices = output / DDM_MATRICES_DIRECTORY
        if destination_matrices.is_symlink() or (
            destination_matrices.exists() and not destination_matrices.is_dir()
        ):
            raise FileExistsError(f"matrix output path is not a directory: {destination_matrices}")
        sources = (video, staging / METADATA_NAME, *matrix_paths)
        destinations = tuple(output / source.relative_to(staging) for source in sources)
        conflicting = next((path for path in destinations if path.is_dir()), None)
        if conflicting is not None:
            raise FileExistsError(f"demo output file path is a directory: {conflicting}")
        destination_matrices.mkdir(exist_ok=True)
        for source, destination in zip(sources, destinations):
            source.replace(destination)
    return output


def main(argv: Sequence[str] | None = None) -> int:
    """Run the demonstration dataset generator."""
    parser = argparse.ArgumentParser(
        prog="ddmsoft-demo",
        description="Generate a deterministic AVI, metadata, and DDM matrix set.",
    )
    parser.add_argument("output", nargs="?", type=Path, help="directory to create or populate")
    parser.add_argument("-o", "--output", dest="output_option", type=Path, help=argparse.SUPPRESS)
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="replace demo-owned files in a nonempty directory without removing other files",
    )
    arguments = parser.parse_args(argv)
    if arguments.output is not None and arguments.output_option is not None:
        parser.error("output must be provided either positionally or with --output, not both")
    requested_output = arguments.output_option or arguments.output
    if requested_output is None:
        parser.error("an output directory is required")
    try:
        output = generate_demo(requested_output, overwrite=arguments.overwrite)
    except (OSError, RuntimeError, ValueError) as error:
        parser.error(str(error))
    print(output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
