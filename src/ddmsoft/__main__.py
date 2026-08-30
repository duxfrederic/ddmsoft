"""Command-line entry point for DDMSoft."""

from __future__ import annotations

import argparse
from collections.abc import Sequence


def main(argv: Sequence[str] | None = None) -> int:
    """Start the native DDMSoft application shell."""
    parser = argparse.ArgumentParser(prog="ddmsoft", description="DDMSoft desktop application")
    parser.add_argument("--version", action="version", version="ddmsoft 0.1.0")
    parser.parse_args(argv)

    from PySide6.QtWidgets import QApplication

    from .gui import create_main_window

    application = QApplication.instance() or QApplication([])
    window = create_main_window()
    window.show()
    return application.exec()


if __name__ == "__main__":
    raise SystemExit(main())
