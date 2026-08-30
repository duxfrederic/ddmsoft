"""Import-safe package shell for the DDMSoft modernization."""

from .models import DDMData, FitRange, FitRequest, FitResult, VideoMetadata
from .io import load_directory, load_matrices, parse_metadata_file, save_matrix_csv
from .engine import compute_ddm, compute_video_ddm
from .fitting import MODEL_REGISTRY, fit_ddm
from .combining import average_ddm, average_groups, merge_ddm

__version__ = "0.1.0"

__all__ = [
    "DDMData",
    "compute_ddm",
    "compute_video_ddm",
    "MODEL_REGISTRY",
    "fit_ddm",
    "average_ddm",
    "average_groups",
    "merge_ddm",
    "FitRange",
    "FitRequest",
    "FitResult",
    "VideoMetadata",
    "load_directory",
    "load_matrices",
    "parse_metadata_file",
    "save_matrix_csv",
    "__version__",
]
