# DDMSoft Legacy Parity Checklist

This inventory describes the committed 2019 application at `2ba8ec6`. It is a
behavior reference for the modernization, not a design proposal. The source
checkout contains only the legacy flat modules and the original `logo/` assets;
no implementation files from the `fix_up` branch are used.

## Baseline

- Starting source: `main` at `2ba8ec6` (`Update windows installation link in README`).
- Working branch at inventory time: `main`.
- Local reference runtime: Python 3.11.13. The modernization CI target is
  Ubuntu 22.04 or newer and Windows 11, each with Python 3.12 and 3.14.
- There is no CI workflow in the baseline repository. Package 1 adds the
  four-entry CI matrix; this document records the baseline rather than claiming
  that those jobs have already run.
- The current checkout has no tracked build or egg metadata. Generated
  `build/`, `ddmsoft.egg-info/`, and `tests/` artifacts mentioned by the
  roadmap remain outside the modernization source and must not be committed.

## Reachable Menus and Actions

These entries are present in the menu declared by `layoutFile.py` and handled by
the event loop in `DDMSoft.py`.

| Menu | Action | Handler / effect |
| --- | --- | --- |
| File | Open directory | Loads AVI files, acquisition parameters, and existing matrices; selects the first matrix when available. Keyboard open shortcuts are also accepted. |
| File | Exit | Stops the main event loop. |
| Tools | Concatenate videos | Selects AVI files and invokes the optional ffmpeg concatenation helper. |
| Tools | Split a video in N matrices (time dependant DDM) | Selects videos, partitions each video, computes one matrix per partition, and refreshes the matrix list. |
| Plotting | Show the DDM matrix (image) | Opens measured and, if present, fitted DDM matrix plots. |
| Plotting | Plot some autocorrelation functions | Plots sampled measured correlation functions and optional fits over the selected q range. |
| Batch, export | Save the DDM matrix as a text file | Writes three tab-delimited CSV files for matrix, lag times, and q values. |
| Batch, export | Save the current fit parameters | Writes a tab-delimited fit text/CSV file, optionally including hydrodynamic radius. |
| Batch, export | Fit and save all the matrices | Fits every discovered matrix and writes one result file per matrix. |
| Batch, export | Save the correlation functions | Writes the refined autocorrelation matrix plus q and lag-time files. |
| More Fitting | CONTIN | Opens the CONTIN parameter window, runs the alpha scan, and opens its plot. |
| Help | About... | Shows the credits dialog. |

## Main Window Controls

The main window preserves this vertical order: computation, fitting, plotting,
then status and progress.

### DDM matrix computation

- Directory input and folder browser (`inputpath`, `browser`), `Load`, and
  `Quit`.
- `Keep existing matrices` / `Re-compute, overwrite existing` radio choices.
- Maximum frame couples (`maxcouples`, default `300`).
- Lag-times-per-decade (`ptperdecade`, default `20`).
- Editable multiline video description containing name, frame rate, and pixel
  size.
- `Process`.
- Direction count (`Nangle`, default `1`). A value above one requests directional
  averaging.

### DDM matrix fitting

- `Merge matrices` and `Average matrices` selection actions.
- Computed matrix selector (`computedlist`).
- Fit model selector (`fitmodel`) with eight legacy models: single exponential,
  second- and third-order cumulants, stretched exponential, double exponential
  with second stretch, exponential with flow, stretched exponential with flow,
  and double exponential with flow.
- q-min and q-max sliders (`qminslider`, `qmaxslider`).
- t-min and t-max sliders (`dtsminslider`, `dtsmaxslider`).
- `Initial guess for the fit`, `Fit the selected matrix`, and `Show the fitted
  parameters`.
- Slider changes force the lower endpoint below the upper endpoint. The legacy
  backend treats upper values as exclusive slice stops; this mismatch is
  explicitly corrected by the roadmap's inclusive range contract in Package 5.

### Plotting

- Optional temperature input in degrees Celsius.
- Optional viscosity input, with the literal `water` selecting the water
  viscosity formula.
- `Plot the matrix and the fit`.
- `Plot the amplitude, the noise, the diffusion`.

## Secondary Dialogs and Interactive Plots

- Initial-guess dialog: one value input and one fixed checkbox per model
  parameter, followed by `submit`.
- CONTIN dialog: q index, minimum/maximum decay rate and count, minimum/maximum
  alpha and count, maximum iterations, `Estimate`, `Plot`, `Save`, and a
  `save all the fits?` checkbox.
- Matrix-selection dialog used by merge and average: text inclusion/exclusion
  filters, one checkbox per matrix, and `Submit`.
- Directory/file dialogs are used for opening directories, saving exports,
  selecting videos, and choosing concatenation output.
- The main DDM plot has a q slider and responds to left/right keys and mouse
  wheel navigation. It updates measured data, fitted data, q text, axis limits,
  and selected time markers.
- The CONTIN plot has an alpha slider and responds to left/right keys and mouse
  wheel navigation. It updates the fit, amplitude/noise, and distribution.
- Matplotlib plots are currently module-global and Tk-backed. This is recorded
  behavior only; Package 8 replaces the implementation with independent Qt
  controllers.

## Dead or Unreachable Legacy Controls

- `Rename the timestamps in a directory` has a handler in `DDMSoft.py` but is
  absent from `menu_def`; it is not reachable from the current menu and is not
  required for parity unless a later review promotes it.
- The POSIX and Windows layout branches duplicate controls with different fixed
  sizes and fonts. Only the branch matching the host OS is constructed.
- The `layoutFile.py` `__main__` block is a layout smoke window, not the
  application launcher; the real application starts at module import and enters
  its event loop.
- The `RadialAverager_test` class is an experimental/reference helper and is not
  called by the active DDM workflow.

## Legacy File Names and Exports

### Matrix persistence

For a video `<directory>/<stem>.avi`, the matrix files are stored in
`<directory>/ddm_matrices/` as:

```text
<stem>_DDM_matrix.npy
<stem>_deltaTs.npy
<stem>_QS.npy
```

The loader discovers a dataset only when all three files exist. Directional
outputs insert an angle token before each suffix, for example
`<stem>_0.0__DDM_matrix.npy`, with one set per sector. Time-dependent outputs
insert a partition token such as `<stem>__i=0__` before the suffix. Existing
matrix names must remain loadable.

Merge output uses `Merged_<title>` plus each of the three suffixes. Average
output uses `Averaged(<frame-rate>)_<title>` plus each suffix. The output
directory is derived from the last selected input matrix.

### Text exports

- Matrix export appends `_DDM_matrix.csv`, `_deltaTs.csv`, and `_QS.csv` to the
  selected base path. All are tab-delimited NumPy text files with `%.6e` values.
- Correlation export writes `<base>_autocorrelationmatrix.csv`, `<base>_qs.csv`,
  and `<base>_dts.csv`, tab-delimited with `%.6e` values. Its matrix contains
  `1 - (DDM - B) / A` for the selected q range.
- Fit export defaults to `.txt` unless the requested path already ends in
  `.txt` or `.csv`. Its first row is tab-delimited and starts with
  `q [m^-1]`, `A`, `B`, followed by model parameter names. Values use `%.3e`.
  When temperature and viscosity are supplied, `R_H eff (nm)` is appended.
- CONTIN export is human-readable tab-delimited text containing matrix, q,
  selected alpha, every candidate alpha, its residual, distribution, and fitted
  correlation data. The current implementation labels the selected candidate
  as the alpha with the smallest residual.

## Visual Reference

The project page's annotated interface image is:

`https://duxfrederic.github.io/ddmsoft/figures/interfaceannotated.png`

The page that embeds it is:

`https://duxfrederic.github.io/ddmsoft/`

These URLs were checked on 2026-08-30. The image is a visual reference only;
the checklist above is the authoritative control and workflow inventory.
