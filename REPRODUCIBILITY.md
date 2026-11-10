# Reproducibility record

This document defines the supported reproduction paths and the evidence boundary for the 2DVSMSI phantom analysis of the accompanying paper. The reference outputs in `results/` (`statistics_results_final.txt` and `numbers*.tex`) are the numerical source of record.

## Scope of the evidence

The analysis concerns one agar phantom containing one titanium rod and one continuous cooling series. PFA, fully sampled 2DVSMSI and GRE/PRF acquisitions were collected at separate scan times during that series. The session acquired nine time points (acquisitions 1-9, numbered in acquisition order). Calibration and LOOCV use all nine, with deviations referenced to each voxel's mean over all nine (a constant reference is absorbed by the ROI intercepts, so slopes and LOOCV do not depend on it). Temperature-change maps and Bland-Altman and Deming analyses use time points 3-9, with references averaged over 3-9. Time point 1 is excluded from that comparison because its GRE scan preceded its 2DVSMSI scans by about 5 min (302 s) instead of 33-36 s; time point 2 because its GRE phase deviates from the cooling trend. Single-time-point illustrations, SNR and display thresholds use time point 3; the GRE magnitude image (series 0016) is from time point 2 because no other GRE magnitude image was stored.

GRE PRF conversion: the phase series is the scanner's phase-difference image of the dual-echo `gre_field_mapping` scan (echoes at 5.00 and 7.46 ms, DICOM EchoNumber 2), so the effective echo-time difference is 2.46 ms. It is defined once in `code/gre_params.m` and used by all PRF code. B0 is the DICOM ImagingFrequency (123.251054 MHz) divided by the proton gamma-bar (42.577478 MHz/T), i.e. 2.8947 T, not the nominal 3 T.

The probe arrays used by the analysis are the mean of the 1 Hz fluoroptic log (`data/probe_log_20250425.csv`) over each scan's data-acquisition window (PFA/Full: first to last ADC time stamp, dummy scans excluded; GRE: DICOM AcquisitionTime + 25.4 s, the DICOM SliceMeasurementDuration; the protocol's total scan time of 27 s includes preparation and is not the data window), with the measured 51 s thermometer-versus-scanner clock offset corrected, computed by `code/build_probe_temperatures.m`. Probe temperatures span 36.8-44.9 C over about 65 min.

Four ROIs were positioned adjacent to the fluoroptic probes on the native grid of each modality. Calibration and LOOCV use these observations from the same series. They quantify internal calibration behaviour and held-out-time-point prediction within that series; they are not independent external validation. No repeat heating experiment, tissue cohort, in-vivo experiment or external test set is included.

## Software environment

The paper analysis was run with MATLAB R2023a, using the Statistics and Machine Learning Toolbox and the Image Processing Toolbox. Python 3.12 (numpy, scipy, h5py; `tests/requirements.txt`) is used for the independent verification. The Pulseq MATLAB toolbox is needed only for `sequence/` and `rev_readout_fov.m` and is not bundled. Record `version` and `ver` with a reproduction attempt; binary-identical output across environments is not assumed.

## Directory configuration

`code/project_paths.m` is the only place where paths are resolved.

- `p.root` is the repository root; `p.code` is `code/`; `p.data` is `data/`; `p.mapvbvd` is `third_party/mapVBVD/`.
- `VSMSI_RAW_DIR` optionally points to the folder with the raw Siemens `.dat` files; otherwise `data/` is used. The GRE DICOM files are always read from `data/`.
- Regenerated files go to `p.output` (`output/`, git-ignored, or `VSMSI_OUTPUT_DIR`): figures in `output/figures`, numbers and statistics in `output/results`, working files of the supporting analyses in `output/supporting` (or `VSMSI_REV_DIR`).
- `p.reference` is `results/`, the committed reference outputs; the code never writes there.

Local absolute paths belong in environment variables, not in `.m` files.

## Input provenance

### Raw reconstruction inputs (not distributed)

`code/reconstruct_all.m` lists the expected file names and the cooling order: PFA 2DVSMSI raw `.dat` files (MID00570, 00576, ..., 00618), fully sampled files (MID00572, 00578, ..., 00620), the GRE phase DICOM files (series 0015, 0017, ..., 0031) and the GRE magnitude DICOM file (series 0016). The raw `.dat` files (about 1.5 GB) are available from the corresponding author on reasonable request. The GRE DICOM files are in `data/` in de-identified form (patient, institution and device fields blanked or replaced, Siemens private CSA headers removed, UIDs regenerated; pixel data and the tags used by the code unchanged).

### Derived analysis inputs (in `data/`)

- `T2MAP_and_img4i_native100.mat`: native-grid 2DVSMSI maps and magnitude data;
- `phaseMAPr_and_magmap.mat`: GRE phase and magnitude data;
- `ROI_and_tempr_resolved.mat`: resolved modality-specific ROI coordinates, probe-temperature arrays (`temp1`/`temp3`/`temprp1`) and oil-reference ROIs;
- `probe_log_20250425.csv`: the fluoroptic log read by `code/build_probe_temperatures.m`;
- `YX.mat`: rod geometry and oil-region display mask;
- `noise_dependence_summary.mat`: pseudo-replica SNR (`snr_replica`, `n_replicas`), written by `code/supporting_analyses/rev_write_dependence_summary.m`.

`ROI_and_tempr_resolved.mat` is the authoritative ROI/temperature input for `statistics_analysis_final.m`.

## Supported pipeline

### A. Full rebuild (needs the raw `.dat` files)

```matlab
addpath code
run code/reconstruct_all.m                 % rewrites data/T2MAP_and_img4i_native100.mat and data/phaseMAPr_and_magmap.mat
run code/statistics_analysis_final.m
run code/paper_figure.m
```

```text
raw 2DVSMSI k-space
  -> mapVBVD -> dataprocessing / dataSOS -> POCS (PFA only, pocsnew.m)
  -> RDNLS apparent-T2 fitting -> T2MAP_and_img4i_native100.mat

GRE DICOM -> phaseMAPr_and_magmap.mat

probe_log_20250425.csv (51 s clock offset corrected) + scan time stamps
  -> build_probe_temperatures.m -> ROI_and_tempr_resolved.mat (temp1, temp3, temprp1)

derived maps + ROI_and_tempr_resolved.mat + YX.mat + noise_dependence_summary.mat
  -> statistics_analysis_final.m -> output/figures/*.png, output/results/statistics_results_final.txt/.mat, numbers.tex
```

### B. Analysis from the derived data (no raw data)

```matlab
addpath code
run code/statistics_analysis_final.m
run code/supporting_analyses/rev_gre_timepoint_residual.m
run code/supporting_analyses/rev_gre_calibration_sensitivity.m
run code/paper_figure.m
```

```bash
py -3.12 -I tests/verify_statistics.py
```

This checks the statistical calculations and figure generation; it does not test raw-data decoding or reconstruction. Run the MATLAB scripts sequentially. The regenerated `output/results/numbers*.tex` and `statistics_results_final.txt` are compared with `results/`. `numbers_snr.tex` is produced only from the raw-data pseudo-replica run (`rev_pseudoreplica.m`, `rev_pseudoreplica_analysis.m`, `rev_write_dependence_summary.m`); its reference copy is in `results/`.

### C. Figure and table map

See the table in `README.md`. Every result number in the paper is a LaTeX macro generated by the code: `statistics_analysis_final.m` writes `numbers.tex`, and the supporting scripts write `numbers_snr.tex`, `numbers_gre_residual.tex` and `numbers_gre_sensitivity.tex` (included by `numbers.tex`).

## POCS validation

`tests/test_pocsnew_reconstruction.m` and `tests/test_pocs_symmetric_indices.m` use synthetic data only. They check data consistency, symmetric-window indices and the absence of the one-pixel shift/phase ramp of the unmodified POCS implementation. `code/supporting_analyses/rev_verify_recon.m` (raw data needed) checks that the per-channel reconstruction used by the pseudo-replica simulation reproduces the stored native-grid maps. Geometry convention: `C_pfa = C_legacy + [0, 1]`, `C_full = C_legacy` (`code/analysis_geometry.m`); image reconstruction uses the standard centred FFT convention.

## Statistical interpretation

`statistics_analysis_final.m` performs per-ROI calibration, descriptive between-ROI spread (coefficient of variation) of the calibration coefficients, time-point-cluster bootstrap confidence intervals, within-series LOOCV, a paired PFA-Full difference with a cluster-bootstrap interval, and defined cross-modality agreement summaries. No F tests or power analysis are used.

- ROI-specific LOOCV (the primary metric) refits one calibration (own slope and intercept) per ROI while holding out one time point. The secondary common-coefficient LOOCV refits the model that produces the temperature maps (one common slope plus ROI intercepts) in every fold; the 95% intervals of the overall LOOCV RMSE are cluster bootstraps over whole time points (9 clusters, approximate, may be narrow).
- The overall RMSE across the four ROI-specific LOOCV RMSE values is an aggregate. It is neither a maximum error nor the error of a predictor using a common slope (see the common-coefficient LOOCV for that).
- Bootstrap resampling keeps the observations of a time point together across ROIs; it does not create independent phantom repetitions.
- Display interpolation does not turn the native-grid measurements into an independent higher-resolution data set.
- The fixed 10-pixel border used by RDNLS for thresholding is not translation-invariant; this does not affect any ROI.

## Scripts that overwrite derived files

Three scripts write into `data/` in place; keep a copy of `data/` before rerunning them:

- `code/build_probe_temperatures.m` overwrites `data/ROI_and_tempr_resolved.mat` (probe arrays, scan times, window statistics; the hand-recorded arrays are kept as `*_manual`).
- `code/supporting_analyses/rev_write_dependence_summary.m` overwrites `data/noise_dependence_summary.mat` and `numbers_snr.tex`; run it only after the pseudo-replica scripts (raw data needed).
- `code/reconstruct_all.m` overwrites `data/T2MAP_and_img4i_native100.mat` (Part A, needs raw .dat files) and `data/phaseMAPr_and_magmap.mat` (Part B, needs only the phase DICOM files in `data/`; set the environment variable `VSMSI_PART_B_ONLY=1` to skip Part A). Part B uses the fixed DICOM rescale (RescaleSlope 2, RescaleIntercept -4096, phase = value/4096*pi), read and asserted per file.

## Output review checklist

After any reconstruction, ROI, temperature-array or statistics change:

- `output/results/statistics_results_final.txt/.mat` and `numbers*.tex` were regenerated;
- `py -3.12 -I tests/verify_statistics.py` passes (125 independently recomputed + 12 consistency checks);
- figures were regenerated;
- probe temperatures are the 1 Hz log means over each scan's acquisition window with the 51 s offset corrected;
- no claim describes within-series LOOCV as independent validation.

## Access and licensing boundaries

Raw scanner data are restricted and available from the corresponding author on reasonable request. The MIT licence in `LICENSE` covers the authors' own code and the derived files. Third-party components keep their own terms:

- `code/pocsnew.m` derives from Michael Völker's POCS implementation (BSD 2-clause), see `third_party/pocs/LICENSE.txt` and `NOTICE.md`;
- `third_party/mapVBVD/` is the Siemens reader of Philipp Ehses and contributors (MIT), see `third_party/mapVBVD/LICENSE`.
