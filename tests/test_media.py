from __future__ import annotations

import subprocess
from pathlib import Path

import pytest

from ddmsoft.media import FFmpegError, FFmpegNotFoundError, concatenate_videos


def test_concatenate_videos_uses_temp_output_and_escaped_concat_paths(tmp_path, monkeypatch):
    first = tmp_path / "first clip's.avi"
    second = tmp_path / "second clip.avi"
    first.touch()
    second.touch()
    calls = []
    concat_contents = []

    class FakeProcess:
        returncode = 0

        def __init__(self, command, **kwargs):
            calls.append((command, kwargs))
            concat_contents.append(Path(command[7]).read_text(encoding="utf-8"))
            Path(command[-1]).write_bytes(b"concatenated")

        def communicate(self, timeout=None):
            return "", ""

    monkeypatch.setattr("ddmsoft.media.subprocess.Popen", FakeProcess)
    output = concatenate_videos((first, second), tmp_path / "combined.avi")

    assert output.read_bytes() == b"concatenated"
    command, kwargs = calls[0]
    assert command[:8] == ["ffmpeg", "-n", "-safe", "0", "-f", "concat", "-i", command[7]]
    assert kwargs == {
        "stdout": subprocess.PIPE,
        "stderr": subprocess.PIPE,
        "text": True,
    }
    assert "file '" in concat_contents[0]
    assert "first clip'\\''s.avi" in concat_contents[0]
    assert not Path(command[7]).exists()


def test_concatenate_videos_reports_missing_ffmpeg(tmp_path, monkeypatch):
    first = tmp_path / "first.avi"
    second = tmp_path / "second.avi"
    first.touch()
    second.touch()

    def missing(*args, **kwargs):
        raise FileNotFoundError("ffmpeg")

    monkeypatch.setattr("ddmsoft.media.subprocess.Popen", missing)
    with pytest.raises(FFmpegNotFoundError, match="not found"):
        concatenate_videos((first, second), tmp_path / "combined.avi")


def test_concatenate_videos_reports_nonzero_ffmpeg_exit(tmp_path, monkeypatch):
    first = tmp_path / "first.avi"
    second = tmp_path / "second.avi"
    first.touch()
    second.touch()

    class FailedProcess:
        returncode = 1

        def __init__(self, *args, **kwargs):
            pass

        def communicate(self, timeout=None):
            return "", "codec failure"

    monkeypatch.setattr("ddmsoft.media.subprocess.Popen", FailedProcess)

    with pytest.raises(FFmpegError, match="codec failure"):
        concatenate_videos((first, second), tmp_path / "combined.avi")


def test_concatenate_videos_terminates_ffmpeg_on_cancellation(tmp_path, monkeypatch):
    first = tmp_path / "first.avi"
    second = tmp_path / "second.avi"
    first.touch()
    second.touch()

    class RunningProcess:
        returncode = None
        terminated = False

        def __init__(self, *args, **kwargs):
            pass

        def communicate(self, timeout=None):
            if self.terminated:
                self.returncode = -15
                return "", ""
            raise subprocess.TimeoutExpired("ffmpeg", timeout)

        def terminate(self):
            self.terminated = True

        def kill(self):
            self.terminated = True

    process = RunningProcess()
    monkeypatch.setattr("ddmsoft.media.subprocess.Popen", lambda *args, **kwargs: process)

    from ddmsoft.media import FFmpegCancelled

    with pytest.raises(FFmpegCancelled):
        concatenate_videos(
            (first, second),
            tmp_path / "combined.avi",
            cancel=lambda: True,
        )
    assert process.terminated
    assert not (tmp_path / "combined.avi").exists()


def test_concatenate_videos_does_not_overwrite_destination_created_during_run(
    tmp_path, monkeypatch
):
    first = tmp_path / "first.avi"
    second = tmp_path / "second.avi"
    output = tmp_path / "combined.avi"
    first.touch()
    second.touch()

    class RacingProcess:
        returncode = 0

        def __init__(self, command, **kwargs):
            Path(command[-1]).write_bytes(b"ffmpeg")

        def communicate(self, timeout=None):
            output.write_bytes(b"other process")
            return "", ""

    monkeypatch.setattr("ddmsoft.media.subprocess.Popen", RacingProcess)

    with pytest.raises(FileExistsError):
        concatenate_videos((first, second), output)
    assert output.read_bytes() == b"other process"
