![logo](https://github.com/duxfrederic/ddmsoft/blob/main/logo/logo.png)

DDMSoft was developed during an internship at RWTH Aachen in the summer of 2019. The project was organized by Jerome J. Crassous and funded by a [BioSoft](http://www.ihrs-biosoft.de/ihrs-biosoft/EN/GSP/GSP_node.html) scholarship. For an introduction and examples see the associated github page: https://duxfrederic.github.io/ddmsoft/.

## What is DDMSoft?
DDMSoft converts differential dynamic microscopy (DDM) videos to DDM matrices and provides fitting and export tools for the resulting data. It is designed for use on the computer that captures the videos, allowing quick feedback during a measurement.

Features include:
- fast conversion of videos to DDM matrices
- directional DDM matrices for anisotropic flow
- time-dependent DDM
- matrix inspection with correlation-function plots
- single, stretched, double, cumulant, and flow relaxation models
- CONTIN relaxation-rate distributions
- matrix merging, averaging, and data export

## Installation
DDMSoft is a native Qt desktop application. It does not use Electron, a web browser, or PySimpleGUI. Python 3.10 through 3.14 are supported; Python 3.14 is the reference environment.

Create an environment and install the project from its directory:
```bash
$ git clone https://github.com/duxfrederic/ddmsoft.git
$ cd ddmsoft
$ python -m venv .venv
$ . .venv/bin/activate             # on Windows: .venv\Scripts\activate
$ python -m pip install -e .
$ ddmsoft
```

The video reader uses OpenCV directly, so the old `sk-video` dependency is no longer required. Existing matrices remain compatible: DDMSoft reads and writes the `ddm_matrices/*_DDM_matrix.npy`, `*_deltaTs.npy`, and `*_QS.npy` files used by the 2019 version.

## Workflow
1. Load a directory containing AVI videos and review or edit the metadata in the video catalog.
2. Process videos in the background. Existing matrices can be kept, or overwritten explicitly.
3. Select a matrix in the Matrix inspector to review correlation functions and DDM data.
4. Choose a relaxation model, edit initial guesses if needed, and fit the selected q/time range.
5. Use the embedded plots and export controls for fitted parameters, correlation functions, or matrix data.

Directional DDM, time-dependent DDM, matrix merging/averaging, and CONTIN remain available in the application.
