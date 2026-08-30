"""Safe external media operations used by the optional GUI tools."""

from __future__ import annotations

import os
import subprocess
from collections.abc import Callable, Iterable
from pathlib import Path
from tempfile import TemporaryDirectory


class FFmpegError(RuntimeError):
    """An ffmpeg invocation failed or returned a non-zero exit status."""


class FFmpegNotFoundError(FFmpegError):
    """The configured ffmpeg executable could not be started."""


class FFmpegCancelled(FFmpegError):
    """A cooperative cancellation request stopped ffmpeg before output commit."""


def concatenate_videos(
    videos: Iterable[str | Path],
    output: str | Path,
    *,
    ffmpeg: str = "ffmpeg",
    overwrite: bool = False,
    cancel: Callable[[], bool] | None = None,
) -> Path:
    """Concatenate videos with a safely quoted ffmpeg concat input file."""
    inputs = tuple(Path(video) for video in videos)
    if len(inputs) < 2:
        raise ValueError("at least two videos are required")
    missing = tuple(path for path in inputs if not path.is_file())
    if missing:
        raise FileNotFoundError(", ".join(str(path) for path in missing))
    destination = Path(output)
    if destination.exists() and not overwrite:
        raise FileExistsError(f"output already exists: {destination}")
    destination.parent.mkdir(parents=True, exist_ok=True)
    with TemporaryDirectory(prefix=".ddmsoft-ffmpeg-", dir=destination.parent) as temporary:
        temporary_directory = Path(temporary)
        concat_file = temporary_directory / "inputs.txt"
        temporary_output = temporary_directory / destination.name
        concat_file.write_text(
            "".join(f"file '{_ffmpeg_quote(path.resolve())}'\n" for path in inputs),
            encoding="utf-8",
        )
        command = [
            ffmpeg,
            "-y" if overwrite else "-n",
            "-safe",
            "0",
            "-f",
            "concat",
            "-i",
            str(concat_file),
            "-c",
            "copy",
            str(temporary_output),
        ]
        try:
            process = subprocess.Popen(
                command,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
            )
        except FileNotFoundError as error:
            raise FFmpegNotFoundError(f"ffmpeg executable was not found: {ffmpeg}") from error
        while True:
            try:
                stdout, stderr = process.communicate(timeout=0.05)
                break
            except subprocess.TimeoutExpired:
                if cancel is None or not cancel():
                    continue
                process.terminate()
                try:
                    process.communicate(timeout=2)
                except subprocess.TimeoutExpired:
                    process.kill()
                    process.communicate()
                raise FFmpegCancelled("video concatenation cancelled")
        if process.returncode != 0:
            details = (stderr or stdout or "").strip()
            raise FFmpegError(
                f"ffmpeg failed with exit code {process.returncode}"
                + (f": {details}" if details else "")
            )
        if not temporary_output.is_file():
            raise FFmpegError("ffmpeg completed without producing the output video")
        if cancel is not None and cancel():
            raise FFmpegCancelled("video concatenation cancelled before output commit")
        if overwrite:
            temporary_output.replace(destination)
        else:
            os.link(temporary_output, destination)
    return destination


def _ffmpeg_quote(path: Path) -> str:
    """Quote a POSIX concat-demuxer path, including embedded apostrophes."""
    return str(path).replace("\\", "\\\\").replace("'", "'\\''")


__all__ = ["FFmpegCancelled", "FFmpegError", "FFmpegNotFoundError", "concatenate_videos"]
