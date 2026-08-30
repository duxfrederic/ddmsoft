# Scientific and Compatibility Changes

This document records intentional differences and clarified contracts in the
modern DDMSoft implementation. The legacy three-file NumPy matrix format and the
legacy fitting model equations remain compatible unless noted below.

## Video input and DDM calculation

### OpenCV reader

Video input now uses OpenCV instead of scikit-video. Decoded BGR/BGRA frames are
converted explicitly to grayscale. The reader validates that the video opens,
that at least one frame decodes, and that every decoded frame has the same
shape. It counts frames actually decoded rather than trusting the container's
reported count, and always releases `VideoCapture`, including on failure.

Codec support is therefore determined by the installed OpenCV build and host
operating system. This is an input compatibility change, not a change to the DDM
power-spectrum definition.

### Short-video lag handling

Lag generation was corrected so every valid video produces usable lags:

- Fewer than two frames is rejected.
- Every lag is a unique, increasing integer in `[1, frame_count - 1]`.
- Lag 1 is always included.
- With the default sampling, videos of two through nine frames use every integer
  lag from 1 through `frame_count - 1` instead of producing no lag values.
- Non-positive or non-integer points-per-decade values are rejected.

### Square frames only

The radial-frequency calibration inherited from the laboratory workflow is
defined only for square frames in this release. Non-square frames are rejected
before FFT calculation rather than being accepted with a silently incorrect q
grid. Supporting rectangular frames in the future requires a separately
validated radial-frequency convention.

### Legacy directional convention

Directional DDM preserves the original formula and inversion-symmetric sector
convention. For Fourier row and column frequencies, the angle is calculated as

```text
theta = atan(q_column / q_row) + pi/2
```

and folded over 180 degrees so opposite Fourier directions share a sector.
Sector 0 is centered on the horizontal Fourier axis; sector centers increase
counter-clockwise in image coordinates.

For `N` sectors, sector `i` is centered at `i*pi/N` and has the wrapped interval

```text
((i - 1/2)*pi/N, (i + 1/2)*pi/N]
```

The lower boundary is excluded and the upper boundary is included. The origin,
whose angle is undefined, contributes to every sector as in the legacy code.
Directional output filenames use the sector center, not an edge:
`0.0`, `180/N`, ..., `(N-1)*180/N` degrees, formatted to one decimal place.
Legacy directional filenames remain discoverable.

### Time-dependent partitioning

Time-dependent DDM now uses contiguous, array-split-equivalent partitions. Every
source frame belongs to exactly one partition, including remainder frames that
the legacy implementation discarded. Zero partitions, more partitions than
frames, and partitions too short for DDM are rejected. Output names contain the
partition's actual zero-based starting frame, for example `__i=7__`.

## Fitting

### Inclusive ranges

`q_min`, `q_max`, `time_min`, and `time_max` are inclusive indices throughout
the GUI, fitting API, and plot markers. Selecting the last q or lag position
includes that sample. Conversion to NumPy's exclusive slice stop occurs once by
adding one to each upper endpoint. The GUI continues to require at least two
positions between each lower and upper selection.

Amplitude and background estimates are now aligned to the same selected q
slice. A nonzero q-min no longer reuses estimates starting at q index zero.

### Structured fit failures

Fitting still uses the legacy model equations and unconstrained Nelder-Mead
optimizer. Each q value now returns its own convergence flag and diagnostic
message. An optimizer error or non-convergence at one q does not falsely mark
the complete fit successful and does not discard later q fits. Curves and
parameters associated with a failed q are fallback diagnostic values and must
not be interpreted as a converged result.

## Matrix combination

Merge and average operations no longer mutate input matrices. They validate
matrix width, q grids, increasing lag grids, overlap scaling, and compatible
shapes before producing output. Averaging groups data by compatible lag-time
arrays and returns one result per group instead of conflating different frame
rates.

## CONTIN

The legacy implementation described alpha selection as an L-curve criterion
but selected the candidate with the smallest residual. This release retains
that behavior for scientific parity and names it accurately as **minimum
residual**. It does not claim to implement L-curve curvature selection. Every
candidate alpha, residual, distribution, fitted curve, amplitude, and noise
value is retained for inspection and export.

Gamma and alpha ranges retain the legacy **linear** spacing, including both
endpoints. They were not silently changed to logarithmic spacing.

## Hydrodynamic conversion

Stokes-Einstein conversions use explicit SI quantities:

- diffusion coefficient: `m^2/s`
- absolute temperature: `K`
- dynamic viscosity: `Pa s`
- internal hydrodynamic radius: `m`

The GUI converts its Celsius temperature field to kelvin and exports radius in
nanometres where labelled. Water viscosity is calculated in `Pa s`. CONTIN
gamma-to-radius conversion uses the same SI contract.

## Excluded unreachable behavior

The legacy source contains a timestamp-renaming handler, but the action was not
present in the reachable menu. It is intentionally excluded from this release
rather than being presented as supported functionality. Existing video and
matrix filenames are not renamed by DDMSoft.
