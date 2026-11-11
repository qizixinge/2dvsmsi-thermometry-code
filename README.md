# 2DVSMSI thermometry near a titanium rod: analysis code and derived data

Code and derived data that reproduce the numbers and figures of

> S. Li, F. Wang, S. Gao, "Sub-minute T2-sensitive multispectral MR thermometry near a metallic implant using view sharing: A phantom feasibility study," *IEEE Transactions on Biomedical Engineering* (citation placeholder: volume, pages, year and DOI to be added after publication).

The study is a phantom proof of concept: one agar phantom with one titanium rod and one continuous cooling series, nine time points. Apparent-T2 maps from partial-Fourier (PFA, POCS) and fully sampled two-dimensional view-shared multispectral imaging (2DVSMSI) are calibrated against fluoroptic probe temperatures and compared with gradient-echo (GRE) proton-resonance-frequency (PRF) thermometry. Calibration and leave-one-out cross-validation (LOOCV) are within the same series; they are not independent external validation.

This repository contains only what is needed to reproduce the reported results: no manuscript sources and no figures. Figures are regenerated into a git-ignored folder.

## Raw data are not included

The raw Siemens k-space files (`meas_MID*.dat`, about 1.5 GB) are **not** in this repository. As stated in the paper, they are available from the corresponding author on reasonable request. Everything downstream of the reconstruction runs from the derived data in `data/`:

| Step | Needs raw `.dat`? |
|---|---|
| Statistics, calibration, LOOCV, bootstrap, Bland-Altman, temperature maps, Figs. 4-5 and S1-S4 (`statistics_analysis_final.m`) | no |
| Image panels of Fig. 3 and Fig. S5 (`paper_figure.m`) | no |
| GRE supporting analyses (`rev_gre_timepoint_residual.m`, `rev_gre_calibration_sensitivity.m`) | no |
| Probe temperatures from the fluoroptic log (`build_probe_temperatures.m`) | no |
| Python verification (`tests/verify_statistics.py`) and MATLAB POCS tests | no |
| GRE phase import (`reconstruct_all.m`, Part B) | no (uses the de-identified DICOM files in `data/`) |
| Reconstruction and T2 fitting (`reconstruct_all.m`, Part A), pseudo-replica SNR (`rev_pseudoreplica*.m`, `rev_write_dependence_summary.m`), `rev_verify_recon.m` | **yes** |

Put raw files into `data/` or point `VSMSI_RAW_DIR` to their folder.

## Requirements

- MATLAB R2023a (used for the paper) with the Statistics and Machine Learning Toolbox and the Image Processing Toolbox. Other releases are untested; record the output of `ver` with a reproduction attempt.
- Python 3.12 with `pip install -r tests/requirements.txt` (numpy, scipy, h5py) for the independent verification.
- Pulseq MATLAB toolbox (https://github.com/pulseq/pulseq), only to run the scripts in `sequence/` or `code/supporting_analyses/rev_readout_fov.m`. It is not bundled, and the exact version used for the scanned sequences was not recorded. The scripts use `mr.opts` with `GradUnit`, `mr.makeExtendedTrapezoid`, `mr.makeAdc` and `seq.write`, so a Pulseq 1.3 or later release is expected.

## Quick start

From the repository root, in MATLAB, one script at a time:

```matlab
addpath code
run code/statistics_analysis_final.m                            % figures, statistics_results_final.txt/.mat, numbers.tex
run code/supporting_analyses/rev_gre_timepoint_residual.m       % numbers_gre_residual.tex
run code/supporting_analyses/rev_gre_calibration_sensitivity.m  % numbers_gre_sensitivity.tex
run code/paper_figure.m                                         % Fig. 3 and Fig. S5 panels
```

then in a shell:

```bash
py -3.12 -I tests/verify_statistics.py     # independent check of the MATLAB results
```

It prints "PASS: 125 independently recomputed checks + 12 consistency checks" (current counts): the first group recomputes values from the inputs (probe temperatures of all 27 acquisition windows from the 1 Hz log, per-ROI and common-coefficient LOOCV, slopes, R², bootstrap intervals, in-sample SDs, B0 and echo times from the DICOM header); the second compares a value stored in one file with its copy in another (for example the LaTeX macros against the .mat file).

Three scripts overwrite derived files in `data/` in place (`code/build_probe_temperatures.m`, `code/supporting_analyses/rev_write_dependence_summary.m`, `code/reconstruct_all.m`); see REPRODUCIBILITY.md, section "Scripts that overwrite derived files", and keep a copy of `data/` before rerunning them.

MATLAB tests: `runtests('tests')`.

### Output convention

The code never overwrites `results/`. Regenerated files go to `output/` (git-ignored): `output/figures/`, `output/results/` and `output/supporting/`. `results/` holds the reference outputs of the paper (`numbers*.tex`, `statistics_results_final.txt/.mat`, `gre_*.mat`, `statistics_independent_verification.json`); compare your `output/results/` with it. `VSMSI_OUTPUT_DIR` moves the output folder. The GRE supporting scripts and `verify_statistics.py` read `statistics_results_final.mat` from `output/results/`, so run `statistics_analysis_final.m` first.

## Reproducing the paper

Every result number in the paper is a LaTeX macro written by the code (`numbers.tex`, which includes `numbers_snr.tex`, `numbers_gre_residual.tex` and `numbers_gre_sensitivity.tex`).

| Paper item | Produced by | Output |
|---|---|---|
| Fig. 3 (A-C), Fig. S5 | `code/paper_figure.m` | `output/figures/` `GRE`, `BVSMSI`, `AVSMSI`, `MIP`, `SOS` (.png and .pdf) |
| Fig. 4, Fig. S1 | `statistics_analysis_final.m` | `Fig4_calibration_SOS.png`, `FigS1_calibration_MIP.png` |
| Fig. 5, Fig. S2 | `statistics_analysis_final.m` | `Fig5_dT_maps.png`, `FigS2_absolute_T_maps.png` |
| Fig. S3 | `statistics_analysis_final.m` | `FigS3_bland_altman.png` |
| Fig. S4 | `statistics_analysis_final.m` | `prfmap.png`, `dt2mapB.png`, `dt2mapA.png` |
| Table 1, supplementary tables, statistics in the text | `statistics_analysis_final.m` | `output/results/statistics_results_final.txt`, `numbers.tex` |
| GRE time-point residual, GRE sensitivity | `rev_gre_timepoint_residual.m`, `rev_gre_calibration_sensitivity.m` | `numbers_gre_residual.tex`, `numbers_gre_sensitivity.tex` |
| SNR (pseudo-replica) | `rev_pseudoreplica.m`, `rev_pseudoreplica_analysis.m`, `rev_write_dependence_summary.m` (raw data needed) | `numbers_snr.tex`, `data/noise_dependence_summary.mat` |
| Readout field of view (265.8 mm) | `rev_readout_fov.m` (Pulseq; `.seq` files written by `sequence/`) | console output, `output/supporting/readout_fov.mat` |
| Fig. 1 (photograph), Fig. 2 (drawing) | not code outputs | not included |

Analysis conventions in brief: nine time points; calibration and LOOCV use all nine; temperature-change maps and Bland-Altman analyses use time points 3-9; uncertainty is descriptive (coefficient of variation) with cluster-bootstrap confidence intervals, including intervals of the overall LOOCV RMSE; the primary LOOCV calibrates each ROI separately and a secondary common-coefficient LOOCV (the model used for the maps) is also reported; the PFA-Full difference is a paired comparison; B0 for the PRF conversion comes from the DICOM header (ImagingFrequency 123.251054 MHz); the GRE echo-time difference is 2.46 ms (`code/gre_params.m`). See `REPRODUCIBILITY.md`.

## Repository layout

| Path | Content |
|---|---|
| `code/` | Reconstruction (`reconstruct_all.m`, `MAP.m`, `dataprocessing.m`, `dataSOS.m`, `pocsnew.m`, `RDNLS.m`), probe temperatures (`build_probe_temperatures.m`), PRF (`PRFreconnew.m`, `gre_params.m`), calibration (`tempfitT2old.m`), `analysis_geometry.m`, `statistics_analysis_final.m`, `paper_figure.m`, `project_paths.m` (single place for all paths) |
| `code/supporting_analyses/` | Scripts that produce reported numbers: GRE residual and sensitivity, pseudo-replica SNR, reconstruction check, readout FOV |
| `sequence/` | Pulseq writers of the two 2DVSMSI protocols (`writemy4TE1b_100_SHARE_full.m`, `write4TE1b_100_SHARE_pocs.m`) |
| `data/` | Derived inputs and GRE DICOM files (see below) |
| `results/` | Reference outputs of the paper |
| `tests/` | `verify_statistics.py`, MATLAB POCS tests, `requirements.txt` |
| `third_party/` | `mapVBVD/` (MIT) and `pocs/` (BSD 2-clause original and a notice on the modified `pocsnew.m`) |

## Data

| File | Content |
|---|---|
| `T2MAP_and_img4i_native100.mat` | native 100x100 2DVSMSI magnitude images and apparent-T2 maps (PFA and Full, SOS and MIP, nine time points) |
| `phaseMAPr_and_magmap.mat` | GRE phase (nine time points, rad) and the single stored GRE magnitude image (time point 2) |
| `ROI_and_tempr_resolved.mat` | ROI coordinates, oil-reference ROIs and probe-temperature arrays (`temp1` PFA, `temp3` Full, `temprp1` GRE) |
| `YX.mat` | rod geometry and oil-region display mask |
| `noise_dependence_summary.mat` | pseudo-replica SNR (`snr_replica`, `n_replicas`) |
| `probe_log_20250425.csv` | 1 Hz fluoroptic log (four channels) of the session |
| `GRE_series00NN_phase.IMA`, `GRE_series0016_magnitude.IMA` | GRE DICOM files (series 0015, 0017, ..., 0031 = acquisitions 1-9, phase-difference images; series 0016 = magnitude, echo 1). De-identified (details below the table). Pixel data, echo times, imaging frequency and acquisition times are unchanged. |


De-identification of the GRE DICOM files: Patient's Name and ID were replaced by PHANTOM and PHANTOM01; birth date, sex, age, size, weight, institution name and address, department, station name, referring and performing physician, accession number and device serial number were emptied; study, procedure-step and series/protocol descriptions were replaced by neutral text (PHANTOM_PULSEQ, gre_field_mapping); the two Siemens CSA header blobs (0029,1010) and (0029,1020) were removed. Other Siemens private elements (acquisition parameters, matrix/FoV/coil text and two private copies of the study date) were kept, and so were all acquisition dates and times, because the code and the paper use them to synchronise the GRE scans with the probe log; this is a phantom scan, but the files are not claimed to meet a formal anonymisation profile. Class-type UIDs (SOP Class, Referenced SOP Class, Media Storage SOP Class, Transfer Syntax, Implementation Class) are unchanged; instance and identity UIDs (study, series, SOP instance, frame of reference, referenced SOP instances) were regenerated with one fixed mapping, so references inside and between files stay consistent. Pixel data, echo times, imaging frequency, pixel bandwidth and rescale values are identical to the originals.

## Third-party code

`third_party/mapVBVD/` (Philipp Ehses and contributors, MIT licence, https://github.com/pehses/mapVBVD) reads Siemens twix files. `code/pocsnew.m` derives from Michael Völker's POCS implementation (BSD 2-clause); see `third_party/pocs/`. The MIT licence in `LICENSE` covers the authors' own code and the released derived data (the `.mat` files, `results/` and the de-identified GRE DICOM files); third-party code keeps its own licences.

## Use of generative AI

Parts of the code were developed with the assistance of generative AI tools (Anthropic Claude via Claude Code; OpenAI ChatGPT/Codex) and were reviewed by the authors.

## Contact

Corresponding author: Song Gao, gaoss@pku.edu.cn.
