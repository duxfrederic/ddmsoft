"""Import-safe package shell for the DDMSoft modernization."""

from .combining import average_ddm, average_groups, merge_ddm
from .contin import CONTINResult, contin_ranges, export_contin, run_contin, run_contin_scan
from .engine import compute_ddm, compute_video_ddm
from .fitting import MODEL_REGISTRY, fit_ddm
from .io import load_directory, load_matrices, parse_metadata_file, save_matrix_csv
from .media import FFmpegCancelled, FFmpegError, FFmpegNotFoundError, concatenate_videos
from .models import DDMData, FitRange, FitRequest, FitResult, VideoMetadata
from .time_dependent import (
    PartitionedDDM,
    compute_time_dependent_ddm,
    compute_time_dependent_videos,
    partition_frame_ranges,
)

__version__ = "0.1.0"

__all__ = [
    "MODEL_REGISTRY",
    "CONTINResult",
    "DDMData",
    "FFmpegCancelled",
    "FFmpegError",
    "FFmpegNotFoundError",
    "FitRange",
    "FitRequest",
    "FitResult",
    "PartitionedDDM",
    "VideoMetadata",
    "__version__",
    "average_ddm",
    "average_groups",
    "compute_ddm",
    "compute_time_dependent_ddm",
    "compute_time_dependent_videos",
    "compute_video_ddm",
    "concatenate_videos",
    "contin_ranges",
    "export_contin",
    "fit_ddm",
    "load_directory",
    "load_matrices",
    "merge_ddm",
    "parse_metadata_file",
    "partition_frame_ranges",
    "run_contin",
    "run_contin_scan",
    "save_matrix_csv",
]
