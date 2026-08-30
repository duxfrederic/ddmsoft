# Modernization Progress

Work packages 0 through 2 are complete. The work was committed separately:

- `3e1b4bc` `docs: inventory legacy behavior`
- `00afc25` `build: add modern package shell`
- `c8dd76c` `test: add backend contracts and fixtures`
- `51ad485` `feat: add metadata and matrix I/O`
- `7d8854f` `feat: add OpenCV DDM engine`

## Package 0

- Added `docs/parity-checklist.md`.
- Documented reachable menus, controls, dialogs, event handlers, dead legacy
  controls, matrix filenames, export formats, CI targets, and the project-page
  interface image URL.
- Legacy source behavior was not changed.

## Package 1

- Added `pyproject.toml` with Python `>=3.12,<3.15` support.
- Added runtime dependencies: NumPy, SciPy, Matplotlib, OpenCV, and PySide6.
- Added development dependencies: pytest, pytest-qt, and Ruff.
- Added the import-safe `src/ddmsoft` package shell and explicit launcher:
  `python -m ddmsoft` and the `ddmsoft` console script.
- Added packaged SVG icon access through `importlib.resources`.
- Added GitHub Actions CI for Ubuntu/Windows and Python 3.12/3.14.
- Legacy flat modules remain in place and are not imported by the launcher.

## Package 2

- Added Qt-free typed records in `src/ddmsoft/models.py`:
  `VideoMetadata`, `DDMData`, `FitRange`, `FitRequest`, and `FitResult`.
- Added explicit-SI water-viscosity and Stokes-Einstein helpers in
  `src/ddmsoft/science.py`.
- Added deterministic test-only frame generators and an independent full-FFT
  DDM reference calculation in `tests/fixtures.py`.
- Added generated matrices for all eight legacy fit model identifiers.
- Added legacy three-file matrix round-trip and incomplete-set fixtures.
- Added 11 tests covering the contracts, fixtures, persistence, and SI
  conversions.

## Package 3

- Added the Qt-free `ddmsoft.io` layer using `pathlib.Path`, UTF-8 context
  managers, and explicit legacy matrix suffixes.
- Added line-oriented acquisition metadata parsing that ignores comments and
  blank lines, splits only on the first colon, and validates `framerate` and
  `pixelsize` independently.
- Preserved common metadata and stem-matched per-video metadata assignment.
  Missing, malformed, ambiguous, incomplete, and invalid matrix data now raise
  descriptive typed I/O exceptions.
- Preserved discovery of legacy three-file matrix sets, including directional
  and partitioned names, and added exact-suffix matrix, correlation, and fit
  exports.
- Added temporary-directory coverage for metadata variants, matrix discovery,
  incomplete sets, and export paths.

## Package 4

- Added `ddmsoft.engine` with an in-memory frame iterator protocol, OpenCV
  decoding, explicit grayscale conversion, frame-shape/nonempty validation,
  and guaranteed `VideoCapture.release()` cleanup.
- Replaced the legacy `scikit-video` input dependency with the OpenCV reader and
  moved Tukey import to `scipy.signal.windows`.
- Added compatibility DDM computation with corrected short-video lag
  generation, progress callbacks, cooperative cancellation, isotropic and
  legacy opposite-direction sector averaging, and computation-before-save
  atomic legacy NumPy output.
- Rectangular frames are explicitly rejected for this milestone; square-frame
  behavior is reference-tested rather than silently using an invalid q grid.
- Added deterministic tests for constant, random, directional, short, invalid,
  rectangular, cancellation, progress, atomic-save, and OpenCV cleanup cases.

## Verification

- `python -m pytest`: 27 passed after packages 3 and 4.
- Package 3 was committed as `51ad485`.
- Package 4 verification passes locally in `7d8854f`.
- Package launcher, packaged resource loading, and Qt-free imports were checked.
- The local environment is Python 3.11.13; the package declares the roadmap's
  Python 3.12 and 3.14 target. Ruff was not installed locally, so CI should run
  the configured Ruff check.

The untracked `roadmap.md` is user-provided and was intentionally not committed.
