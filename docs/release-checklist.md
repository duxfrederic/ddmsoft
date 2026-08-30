# Release Checklist

This checklist separates automated evidence from manual release acceptance.
Passing a cited test does not complete the corresponding human workflow step.
Every item below remains **PENDING** until a reviewer performs it with the
installed release candidate and records platform, Python version, media/data,
date, and result.

## Release qualification

- **COMPLETE locally (Linux, Python 3.14.7, 2026-08-30):** Built the sdist and
  wheel, installed the wheel with all dependencies into a fresh virtual
  environment, imported resources from outside the source tree, started an
  offscreen native window, and exercised both installed console tools.
- **COMPLETE locally:** Full unit, integration, and offscreen GUI suite, Ruff,
  bytecode compilation, import/resource side-effect checks, retired dependency
  audit, and bare-exception review.
- **COMPLETE locally:** Automated coverage for spaces and non-ASCII paths,
  read-only inputs, unwritable outputs, corrupt video/matrices, missing metadata,
  and missing ffmpeg.
- **COMPLETE locally:** Representative measurements are recorded in
  [`benchmark-results.json`](benchmark-results.json).
- **COMPLETE locally:** `ddmsoft-demo` generated a decodable AVI, valid metadata,
  and a discoverable legacy matrix set from the installed wheel.
- **PENDING in CI/manual review:** Windows and Python 3.12 verification, native
  Windows path/permission behavior, real laboratory codecs and matrices, and
  the complete human workflow below. The CI workflow defines Linux/Windows and
  Python 3.12/3.14 jobs for this exact release process.

Suggested release commands:

```bash
python -m pytest
python -m ruff check .
python -m build
ddmsoft-demo
ddmsoft-benchmark
```

## Manual workflow mapping

| # | Manual acceptance step | Automated test evidence | Human signoff |
| ---: | --- | --- | --- |
| 1 | Start DDMSoft from an arbitrary working directory. | `test_import_is_side_effect_free_and_resource_loads_outside_repository`; clean-wheel offscreen launch recorded above. | **PENDING:** visible clean-environment launch on each supported OS. |
| 2 | Open a directory containing multiple videos and metadata files. | `test_metadata_supports_one_file_per_video_and_rejects_ambiguous_files`; `test_load_directory_populates_metadata_and_matrix_catalog`. | **PENDING:** laboratory directory and native file dialog review. |
| 3 | Read and edit frame rate and pixel size in a table. | `test_load_directory_populates_metadata_and_matrix_catalog`; `test_invalid_metadata_edits_are_retained_and_block_processing`. | **PENDING:** visual table editing and validation review. |
| 4 | Keep an existing matrix while processing a new video. | `test_video_job_keeps_existing_sets_and_requires_recompute_to_replace`. | **PENDING:** real-video workflow and resulting catalog review. |
| 5 | Recompute a selected existing matrix after explicit confirmation. | The replace behavior is covered by `test_video_job_keeps_existing_sets_and_requires_recompute_to_replace`; no dedicated confirmation-dialog test. | **PENDING:** confirmation wording, default choice, and overwrite result. |
| 6 | Compute isotropic DDM and observe responsive progress. | `test_isotropic_output_matches_independent_full_fft_reference`; `test_processing_uses_worker_and_refreshes_catalog_without_blocking`; `test_worker_keeps_qt_event_loop_responsive_and_reports_monotonic_progress`. | **PENDING:** real codec, progress visibility, responsiveness, and laboratory result review. |
| 7 | Compute multi-sector directional DDM. | `test_oriented_grating_has_independent_directional_reference_response`; `test_directional_processing_discovers_selects_plots_and_fits_each_sector`. | **PENDING:** visual angle convention and laboratory directional-data review. |
| 8 | Select each resulting matrix from the normal selector. | `test_directional_processing_discovers_selects_plots_and_fits_each_sector`; `test_matrix_selection_clamps_ranges_and_disables_without_selection`. | **PENDING:** selector labels and all generated outputs. |
| 9 | Open the main correlation/DDM plot. | `test_fit_workflow_uses_worker_range_and_keeps_modeless_plots_independent`; `test_matrix_plot_contains_measured_and_fitted_artists`. | **PENDING:** visual axes, labels, scaling, and overlays. |
| 10 | Move the plot q slider with mouse, wheel, and arrow keys. | `test_correlation_plot_q_callback_updates_owned_artists` covers programmatic selection, arrow-key, and wheel callbacks. | **PENDING:** physical slider, mouse, wheel, and keyboard interaction. |
| 11 | Keep the plot open while changing main-window q and time sliders. | `test_fit_workflow_uses_worker_range_and_keeps_modeless_plots_independent`. | **PENDING:** visual/modeless interaction review. |
| 12 | Edit initial fit guesses and fix selected parameters. | `test_initial_guess_dialog_is_generated_from_model_and_validates_values`; `test_fixed_values_and_caller_owned_inputs_are_preserved`. | **PENDING:** dialog usability and full model-by-model review. |
| 13 | Fit the selected inclusive range. | `test_inclusive_final_q_and_time_positions_are_fitted`; `test_fit_job_retains_matrix_identity_and_inclusive_request`; `test_fit_workflow_uses_worker_range_and_keeps_modeless_plots_independent`. | **PENDING:** laboratory fit and endpoint confirmation. |
| 14 | Inspect fit overlays and fitted-parameter plots. | `test_matrix_plot_contains_measured_and_fitted_artists`; `test_fit_workflow_uses_worker_range_and_keeps_modeless_plots_independent`. | **PENDING:** scientific and visual interpretation by a reviewer. |
| 15 | Enter temperature and viscosity and inspect hydrodynamic radius output. | `test_water_viscosity_uses_kelvin_and_pa_seconds`; `test_stokes_einstein_round_trips_in_si_units`; `test_contin_radius_conversion_is_explicit_si`. | **PENDING:** GUI entry, units, labels, and known-sample laboratory check. |
| 16 | Export matrix data, correlation functions, and fit parameters. | `test_exports_write_expected_dimensions_and_headers`; `test_exports_use_exact_suffixes_and_utf8_text`. | **PENDING:** open exported files in the laboratory's downstream tools. |
| 17 | Merge compatible matrices and reject incompatible matrices clearly. | `test_merge_refreshes_matrix_catalog_after_successful_write`; `test_combination_rejects_incompatible_axes_and_zero_scaling`. | **PENDING:** manual compatibility/error-message review. |
| 18 | Average multiple compatible frame-rate groups. | `test_average_groups_by_lag_grid_and_names_results`; `test_average_refreshes_matrix_catalog_and_rejects_incompatible_inputs`. | **PENDING:** multi-group GUI output and scientific review. |
| 19 | Run batch fitting and observe per-matrix progress. | `test_batch_fit_exports_each_matrix_and_prevents_silent_overwrite`; `test_batch_fit_can_continue_after_one_matrix_failure`; `test_batch_fit_main_window_stores_results_by_full_path`. | **PENDING:** visible progress, continue/stop policy, and output review. |
| 20 | Split a video into time-dependent matrices without losing frames. | `test_partition_frame_ranges_cover_every_source_frame_without_empty_ranges`; `test_time_dependent_ddm_preserves_actual_partition_starts_and_reports_progress`; `test_time_dependent_main_window_runs_selected_videos_and_refreshes_catalog`. | **PENDING:** real-video partition and frame-accounting review. |
| 21 | Run CONTIN for one q value. | `test_contin_worker_returns_q_identity_and_scans_all_requested_candidates`; `test_contin_returns_all_candidates_and_minimum_residual_selection`. | **PENDING:** laboratory parameter selection, runtime, and result review. |
| 22 | Explore all CONTIN alpha candidates with the plot slider. | `test_main_plot_windows_and_contin_windows_are_independent`; `test_contin_returns_all_candidates_and_minimum_residual_selection`. | **PENDING:** slider, key/wheel behavior, labels, and candidate plots. |
| 23 | Export the chosen or all CONTIN candidates. | `test_contin_export_uses_each_candidate_amplitude_and_noise`; `test_contin_main_window_keeps_independent_result_and_exports_candidates`. | **PENDING:** manually verify selected-only and all-candidate files. |
| 24 | Open multiple independent plot windows. | `test_main_plot_windows_and_contin_windows_are_independent`; `test_fit_workflow_uses_worker_range_and_keeps_modeless_plots_independent`. | **PENDING:** window focus, lifetime, interaction, and DPI review. |
| 25 | Cancel a long-running operation safely. | `test_worker_cancellation_emits_cancelled_without_result`; `test_processing_cancel_button_requests_cooperative_cancellation`; `test_video_job_cancellation_never_commits_partial_outputs`; `test_concatenate_videos_terminates_ffmpeg_on_cancellation`. | **PENDING:** cancel real DDM, fitting, CONTIN, time-dependent, and ffmpeg jobs. |
| 26 | Close the application with and without an active job. | `test_closing_window_waits_for_worker_cancellation`; normal test teardown closes inactive windows. | **PENDING:** native close behavior, prompts, process exit, and output integrity. |

## Signoff record

Record one row per completed walkthrough. Do not mark the release accepted from
automated evidence alone.

| OS | Python | DDMSoft build | Dataset/video | Reviewer | Date | Result and issues |
| --- | --- | --- | --- | --- | --- | --- |
| _pending_ | _pending_ | _pending_ | _pending_ | _pending_ | _pending_ | _pending_ |
