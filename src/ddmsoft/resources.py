"""Access packaged application resources without relying on the working directory."""

from __future__ import annotations

from importlib.resources import files
from importlib.resources.abc import Traversable


def resource(name: str) -> Traversable:
    """Return a packaged resource such as the application icon."""
    return files("ddmsoft").joinpath("resources", name)


def icon() -> Traversable:
    """Return the packaged application icon."""
    return resource("ddmsoft.svg")
