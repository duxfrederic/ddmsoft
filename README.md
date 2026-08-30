# DDMSoft

DDMSoft is a native desktop application for differential dynamic microscopy
(DDM): it reads microscopy videos, computes DDM matrices, fits relaxation
models, and provides interactive plots and exports.

## Install and launch

DDMSoft supports Windows and Linux with Python 3.12, 3.13, or 3.14. A
standalone installer is not required.

```bash
python -m pip install ddmsoft
ddmsoft
```

To install a source checkout or release archive instead:

```bash
python -m pip install .
ddmsoft
```

Using a virtual environment is recommended. `pip` installs the Qt, NumPy,
SciPy, Matplotlib, and OpenCV dependencies.

## Input directory

Open a directory containing AVI videos and UTF-8 acquisition metadata. Use one
shared `.txt` file for all videos, or one stem-matched file per video. Required
values are frame rate in frames per second and pixel size in metres per pixel:

```text
framerate: 30
pixelsize: 1e-6
```

Existing legacy matrix sets in `ddm_matrices/` load without conversion when all
three files are present: `*_DDM_matrix.npy`, `*_deltaTs.npy`, and `*_QS.npy`.

## Primary workflow

1. Open the acquisition directory and review or edit frame rate and pixel size
   in the video table.
2. Choose whether to keep existing matrices or explicitly recompute them, then
   select lag sampling, frame-pair, and direction settings and click `Process`.
3. Select a matrix and inspect its DDM/correlation plot. Use the q and time
   sliders to choose the fit region.
4. Select a model, edit initial guesses or fixed parameters if needed, and run
   the fit.
5. Inspect matrix, correlation, fitted-parameter, amplitude, noise, diffusion,
   and optional hydrodynamic-radius plots; export the required data.

The q-min, q-max, time-min, and time-max slider values are **inclusive**. The
samples at both selected endpoints participate in fitting and plot markers.

## Advanced workflows

- Set more than one direction part to compute legacy opposite-direction sectors
  over 180 degrees. Each sector is saved and selected as a normal matrix.
- Use `Tools` to split videos into time-dependent matrices. Every source frame
  is assigned to exactly one partition.
- Merge compatible lag ranges, average compatible lag-time groups, or batch-fit
  selected matrices from the `DDM matrix fitting` and `Batch, export` actions.
- Use `More Fitting > CONTIN` to scan one q value, inspect every alpha candidate,
  and export the selected candidate or the complete scan.
- Use `Tools > Concatenate videos` for optional ffmpeg-backed concatenation.

## Video codecs and ffmpeg

Video decoding uses OpenCV. Codec availability depends on the OpenCV build and
operating system, so an `.avi` extension does not guarantee that a stream can be
decoded. DDMSoft reports videos that cannot be opened, have no decodable frames,
or change frame shape. The current DDM calculation accepts square frames only;
crop or transcode rectangular video before processing.

Concatenation requires an external `ffmpeg` executable on `PATH`; it is not a
Python dependency. The operation uses stream copy rather than re-encoding, so
the selected videos must have compatible streams and container parameters.
Missing executables and ffmpeg failures are reported without committing a
partial output.

## Demo and benchmark

Generate a small synthetic demonstration dataset:

```bash
ddmsoft-demo ./ddmsoft-demo-data
```

Run the representative synthetic benchmark:

```bash
ddmsoft-benchmark --output ./ddmsoft-benchmark.json
```

These console scripts are included in the release tooling. The current release
qualification measurements are recorded in
[`docs/benchmark-results.json`](docs/benchmark-results.json).

## Release notes

- [Intentional scientific and compatibility changes](docs/scientific-changes.md)
- [Release qualification and manual checklist](docs/release-checklist.md)
- [Legacy parity inventory](docs/parity-checklist.md)
