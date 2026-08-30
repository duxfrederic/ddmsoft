"""Command-line entry point for DDMSoft."""

from __future__ import annotations

import argparse
from collections.abc import Sequence
from importlib.resources import as_file


def main(argv: Sequence[str] | None = None) -> int:
    """Start the native DDMSoft application shell."""
    parser = argparse.ArgumentParser(prog="ddmsoft", description="DDMSoft desktop application")
    parser.add_argument("--version", action="version", version="ddmsoft 0.1.0")
    parser.parse_args(argv)

    from PySide6.QtGui import QIcon
    from PySide6.QtWidgets import QApplication

    from .gui import create_main_window
    from .resources import icon

    application = QApplication.instance() or QApplication([])
    with as_file(icon()) as icon_path:
        application.setWindowIcon(QIcon(str(icon_path)))
        window = create_main_window()
        window.show()
        return application.exec()


if __name__ == "__main__":
    raise SystemExit(main())
