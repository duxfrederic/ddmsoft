"""Import-safe package shell for the DDMSoft modernization."""

from .models import DDMData, FitRange, FitRequest, FitResult, VideoMetadata

__version__ = "0.1.0"

__all__ = [
    "DDMData",
    "FitRange",
    "FitRequest",
    "FitResult",
    "VideoMetadata",
    "__version__",
]
