# Modernization Progress

Work packages 0 through 7 are complete. The work was committed separately:

- `3e1b4bc` `docs: inventory legacy behavior`
- `00afc25` `build: add modern package shell`
- `c8dd76c` `test: add backend contracts and fixtures`
- `51ad485` `feat: add metadata and matrix I/O`
- `7d8854f` `feat: add OpenCV DDM engine`
- `f0884eb` `feat: add structured fitting API`
- `46108f8` `feat: add pure DDM combination workflows`
- `ad9a5b2` `feat: isolate CONTIN backend`

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

## Package 5

- Added the Qt-free `ddmsoft.fitting` API and a single model registry for all
  eight legacy model identifiers, including stable IDs, display labels,
  functions, parameter definitions, and defaults.
- Added structured fitting from `DDMData` and `FitRequest` with inclusive q and
  time ranges, per-q A/B estimates aligned to the selected q slice, and
  returned correlation/fitted curves generated from the returned parameters.
- Preserved Nelder-Mead fitting and legacy model formulas while scaling optimizer
  coordinates for the existing SI-valued parameters.
- Added per-q convergence status and messages, continuation after an optimizer
  failure, progress/cancellation hooks, and caller-input immutability.
- Corrected the deterministic cumulant fixture parameters to use the SI scale
  implied by their model formula, keeping generated data numerically meaningful.
- Added generated-data coverage for every model, inclusive final positions,
  nonzero q-min alignment, fixed values, failure continuation, and cancellation.

## Package 6

- Added pure `DDMData` merge and average operations with copied output arrays,
  q-grid, matrix-shape, lag-order, overlap, and zero-scaling validation.
- Added deterministic compatible-lag grouping with one named average result per
  group; grouping uses lag arrays rather than inferred frame-rate floats.
- Added array-split-equivalent time-dependent partition ranges that account for
  every frame and retain each partition's actual start/stop frame metadata.
- Added sequential partition computation and per-video metadata dispatch, plus
  I/O-layer saving with legacy partition name tokens.
- Added tests for input immutability, grouping, invalid axes/lags, zero scaling,
  frame accounting, per-video metadata, and partition output names.

## Package 7

- Added the Qt-free `ddmsoft.contin` result API, replacing generator-style
  progress/final-result handling with `run_contin` and a structured
  `CONTINResult` containing every alpha candidate.
- Added validated gamma/alpha min/max/count helpers, progress and cooperative
  cancellation between alpha and optimizer iterations, and minimum-residual
  selection labeled accurately rather than as an L-curve criterion.
- Removed GUI-window coupling from the modern CONTIN exporter and corrected
  candidate-specific amplitude/noise export values.
- Added optional SI diffusion-rate to hydrodynamic-radius conversion and tests
  for candidate retention, validation, cancellation, export, and conversion.

## Verification

- `python -m pytest`: 57 passed after packages 3 through 7.
- Package 3 was committed as `51ad485`.
- Package 4 verification passes locally in `7d8854f`.
- Package 5 verification passes locally in `f0884eb`.
- Package 6 verification passes locally in `46108f8`.
- Package 7 verification passes locally in `ad9a5b2`.
- Package launcher, packaged resource loading, and Qt-free imports were checked.
- The local environment is Python 3.11.13; the package declares the roadmap's
  Python 3.12 and 3.14 target. Ruff was not installed locally, so CI should run
  the configured Ruff check.

The untracked `roadmap.md` is user-provided and was intentionally not committed.
