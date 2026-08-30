"""Command-line entry point for DDMSoft."""

from __future__ import annotations

import argparse
from collections.abc import Sequence


def main(argv: Sequence[str] | None = None) -> int:
    """Run the temporary package launcher without starting a GUI."""
    parser = argparse.ArgumentParser(prog="ddmsoft", description="DDMSoft desktop application")
    parser.parse_args(argv)
    print("DDMSoft is installed; the native application shell is not available yet.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
