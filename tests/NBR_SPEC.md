# ROI neighbourhood analysis — specification (2026-11-27)

Question: how much do the reported calibration slopes, LOOCV RMSEs and between-ROI slope CVs depend on the
exact position of the 3x3 ROI windows? Each ROI window is shifted on its own grid and everything is recomputed.

## Variants, grids, inputs
Five variants, in this order everywhere (index v = 1..5):

| v | name | grid | map | archived 3x3 centres (row, col; 1-based) | probe temperatures (9x4) |
|---|------|------|-----|------------------------------------------|--------------------------|
| 1 | GRE      | 128 | drift-corrected PRF signal (as built in tests/verify_statistics.py for the cfg[0] check) | statistics_results_final.mat `CPRF` | ROI_and_tempr_resolved.mat `temprp1` |
| 2 | PFA-SOS  | 100 | `T2MAPsos(:,:,:,1)` | `C2D` | `temp1` |
| 3 | Full-SOS | 100 | `T2MAPsos(:,:,:,3)` | `C2D_full` | `temp3` |
| 4 | PFA-MIP  | 100 | `T2MAP(:,:,:,1)` | `C2D` | `temp1` |
| 5 | Full-MIP | 100 | `T2MAP(:,:,:,3)` | `C2D_full` | `temp3` |

T2 maps: data/T2MAP_and_img4i_native100.mat (v7.3; MATLAB dims row x col x timepoint x variant, 100x100x9x3).
The ROI position is the same at all nine time points.

## Shifted windows
Offsets (dr, dc) with dr, dc in {-2,-1,0,1,2}. Radius 1 = offsets with |dr|,|dc| <= 1 (9 windows per ROI);
radius 2 = |dr|,|dc| <= 2 (25 windows per ROI). The window for ROI k is rows c_k(1)+dr+(-1:1), cols c_k(2)+dc+(-1:1).

A window is **kept** only by rules that do not use the probe temperatures:
1. all 9 voxels inside the image (1..N);
2. it does not overlap the rod box: overlap means some window row is in [R1,R2] AND some window column is in [C1,C2],
   with [R1 R2 C1 C2] = `geometry.rod_100` (100 grid) or `geometry.rod_128` (128 grid) from statistics_results_final.mat;
3. 2DVSMSI variants only: all 9 voxels finite and > 0 at all 9 time points in that variant's own map
   (SOS map for SOS, MIP map for MIP).

## Per window
- ROI signal y_t (t = 1..9): mean over the 9 voxels of the map at time point t (for GRE: of the PRF signal, same sign
  and scaling as the published GRE calibration, so that the published alpha is reproduced).
- slope: ordinary least squares of y on the ROI's probe temperature (own intercept).
- per-ROI LOOCV RMSE: leave one time point out, refit slope and intercept on the other 8, convert the held-out y to a
  temperature, error = estimate - probe; RMSE over the 9 held-out errors.
- Acceptance: at offset (0,0) the slope and LOOCV RMSE of every variant and ROI equal the archived values in
  statistics_results_final.mat (`ALPHA_PAPER` for GRE in units of 1e-8/degC, `m2_p_sos`, `m2_f_sos`, `m2_p_mip`,
  `m2_f_mip`; `loocv_rmse` rows = variants in the order above) to 1e-9, and all archived windows are kept.

## Combinations (2DVSMSI variants v = 2..5)
A combination = one kept window for each of the four ROIs (Cartesian product of the four kept lists).
For each combination:
- CV = 100 * sample SD (n-1) of the four slopes / |mean of the four slopes| (all means are positive in this data set;
  a combination with mean 0 would have an undefined CV and must be reported, not dropped);
- overall RMSE = sqrt(mean of the four per-ROI LOOCV RMSE^2) (= RMS of the 36 held-out errors);
- common-slope LOOCV RMSE: model y_kt = a_k + m * (T_kt - mean_t T_k); leave one time point out for all four ROIs
  together, refit m and a_k on the other 8, estimate = (y_kt - mean_train y_k)/m + mean_train T_k, error = estimate - T_kt;
  RMS over the 36 held-out errors. (At (0,0) this equals sqrt(mean(common_rmse(v-1,:).^2)) in the results file.)
- Acceptance: at the archived combination, CV and common-slope RMSE equal the published values (1e-9).

Percentile of the archived value x0 within a set x of n values (n includes x0):
pct = 100 * (#{x < x0} + 0.5 * (#{x == x0} - 1)) / (n - 1).
0 = archived is the lowest, 50 = typical.

SOS-versus-MIP share, per acquisition (PFA: v=2 vs 4; Full: v=3 vs 5): over the combinations that are kept in BOTH
variants (same offsets for every ROI), the percentage with CV_SOS < CV_MIP.

Radius 1: all three quantities. Radius 2: CV only (number of combinations is large).

## Outputs (LaTeX macro names without the leading \n; format in brackets)
Radius 1:
- NbrNumWin [%d]: windows kept per ROI (must be 9 for every 2DVSMSI variant and ROI; assert)
- NbrNumComb [%d]: combinations per 2DVSMSI variant (must be equal for all four; assert)
- for V in PfaSos, FullSos, PfaMip, FullMip:
  NbrCvMed{V} [%.1f], NbrCvPct{V} [%.0f], NbrRmsMed{V} [%.2f], NbrRmsPct{V} [%.0f],
  NbrCommonRmsMed{V} [%.2f], NbrCommonRmsPct{V} [%.0f]
  (Med = median over the combinations, Pct = percentile of the archived combination)
- NbrSosLtMipPfa, NbrSosLtMipFull [%.0f]
- GRE per ROI (k = I, II, III, IV): NbrGreRmseMed{k} [%.2f] (median per-ROI LOOCV RMSE over the kept windows; all 9 must
  be kept, assert), NbrGreRank{k} [%d] (1 + number of kept windows with a lower LOOCV RMSE than the archived window)
Radius 2:
- NbrNumCombTwoPfa, NbrNumCombTwoFull [%d]: combinations kept in both SOS and MIP
- NbrCvMedTwo{V} [%.1f] for the four 2DVSMSI variants (median over that variant's own kept combinations)
- NbrSosLtMipTwoPfa, NbrSosLtMipTwoFull [%.0f]

## Shared inputs (not independently verified)
Both implementations take the ROI centres and the rod boxes from statistics_results_final.mat (`geometry`); an error
in these coordinates would be shared. They are checked elsewhere only through the reproduction of the archived results.

## Stored in statistics_results_final.mat (struct NBR, written by statistics_analysis_final.m)
- NBR.off = -2:2
- NBR.slope, NBR.rmse: 5 x 4 x 5 x 5 [variant, ROI, dr index, dc index]; NaN where the window is not kept
  (GRE slope in 1e-8/degC)
- NBR.valid: same size, logical
- NBR.cv, NBR.rms, NBR.crms: 1 x 4 cell (variants 2..5), radius-1 vectors over the combinations (order not specified;
  compare sorted)
