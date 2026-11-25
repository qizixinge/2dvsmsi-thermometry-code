% Calibrate native-grid T2 maps, quantify within-series errors, and export manuscript figures and statistics.
%% statistics_analysis_final.m
% Integrated statistical analysis and figure generation for the T2-based
% 2DVSMSI MR thermometry manuscript (IEEE TBME submission).
%
% WHAT CHANGED IN THIS REVISION (2026-09-03)
%   * All 2DVSMSI quantities are computed on the NATIVE 100x100 grid
%     (T2MAP_and_img4i_native100.mat).  Fit first, interpolate only for
%     display.  See rerun_native_100.m.
%   * Archived ROI coordinates are read from ROI_and_tempr_resolved.mat.
%     analysis_geometry maps them to the corrected centered-FFT frame while
%     preserving the physical sampling regions. Full coordinates are unchanged.
%   * Common-slope calibration for SOS with a cluster bootstrap CI (time points
%     resampled as intact units).  MIP is NOT pooled (see REVISION 2026-11-10).
%   * Per-ROI RMSE / MAE / bias / max|e|, as requested in the 08-28 comments.
%   * Main-text Figure 4 shows PRF + SOS only; MIP curves go to Supplementary
%     Figure S1.  Main-text Figure 5 shows dT maps with invalid voxels as NaN;
%     absolute temperature maps go to Supplementary Figure S2.
%   * Bland-Altman is computed on the native grid and reports bias + LOA only
%     (no proportional-bias p value), and goes to Supplementary Figure S3.
%
% REVISION 2026-10-12 (noise_dependence_summary.mat)
%   * SNR: the upper-left patch lies in the readout roll-off zone and contains
%     a smooth intensity gradient, so its SD is not the image noise SD. The
%     2DVSMSI SNR now uses the pseudo-replica noise SD in the reference
%     rectangle; GRE uses a Rayleigh-corrected background estimate. The legacy
%     patch-SD values are still printed for traceability.
%   * The native voxel is 2.6585 mm (readout) x 2.56 mm (phase encoding)
%     (analysis_geometry.px_100). Rod distances use it, and the native maps are
%     resampled onto the GRE grid by physical position (centred grids) instead
%     of the corner-aligned 1..100 -> 1..128 stretch. This moves the Bland-
%     Altman pairs, Figure 5, Figure S2 and the display oil mask.
%
% REVISION 2026-10-22 (TBME only)
%   * Calibration error bars are +/-1 SD (N-1 normalisation) across the nine
%     voxels of each ROI.
%   * Figures are written to output/figures (project.figures) under the names
%     used by tbme_main.tex and tbme_supplementary.tex; the text and .mat
%     results go to output/results (project.results).
%   * Display masks: PFA and Full each take the threshold from their own
%     time-point-TP_CMP (= 3) 16 ms SOS magnitude (3 x the maximum of the 10-pixel
%     edge band); GRE uses its own magnitude. Previously the PFA value was
%     applied to both 2DVSMSI variants. The masks affect display and the
%     Bland-Altman voxel selection only.%
% REVISION 2026-11-10 (second audit round)
%   * S2b: secondary common-coefficient LOOCV (the model used for the maps), 95% cluster-bootstrap intervals
%     of the overall LOOCV RMSE, in-sample SD of (fit-derived - probe temperature) per ROI (Weber et al. 2017
%     definition) and the cooling rate; all written as macros.
%   * GRE phase uses the fixed DICOM rescale (reconstruct_all.m Part B); GRE probe window 25.4 s.

%
% REVISION 2026-11-05 (nine time points, LaTeX macros)
%   * Time points 1-9 = acquisitions 1-9 of 2025-04-25 (acquisition 1 re-included;
%     every scan has its own probe temperature, see build_probe_temperatures.m).
%   * Calibration, cluster bootstrap and LOOCV use all nine time points.  Deviations are referenced to each voxel's
%     mean over all nine; S1a checks that referencing time points 3-9 instead
%     changes neither the slopes nor the LOOCV errors.
%   * Cross-modality comparisons (dT maps, Bland-Altman, Deming) use time points
%     3-9 (BA_TPS): time point 1 because its GRE preceded its 2DVSMSI scans by
%     ~5 min, time point 2 because its GRE phase image is the known outlier.
%     Single-time-point illustrations (Figs 5, S2, S4, display thresholds, SNR)
%     use TP_CMP = 3.
%   * Every reported number is written as a LaTeX macro to output/results/numbers.tex
%     (S7).  SNR, GRE residual SDs and GRE sensitivity numbers come from
%     numbers_snr.tex / numbers_gre_residual.tex / numbers_gre_sensitivity.tex,
%     written by the supporting scripts and \InputIfFileExists-ed by numbers.tex.
%
% REVISION 2026-11-10 (statistics simplified; single-phantom feasibility study)
%   * The slope-homogeneity F tests, the noncentral-F detectable-slope/power analysis and the
%     "POOL licensed = p > 0.05" gate are REMOVED: they assume independent ROI residuals, which the
%     design does not support (residual correlations up to +0.67 between ROIs).  Replaced by
%     descriptive output only: per-ROI slopes with R^2, between-ROI slope CV, and for SOS the common
%     slope with a cluster-bootstrap CI (whole time points resampled, so ROI independence is not assumed).
%     Descriptive rule (no threshold): a common coefficient is reported for the SOS combination, the
%     primary combination; MIP is reported per ROI because its between-ROI slope CVs were much larger
%     (both printed below).
%   * PFA-vs-Full coefficient difference: paired cluster bootstrap (same resampled time points for
%     both variants, same seed), reported as difference and 95% CI only.
%   * Deming slope of 2DVSMSI dT on GRE dT with the error-variance ratio estimated by Grubbs' three-instrument
%     estimator (PFA, Full, GRE in the same voxels); lambda = 1 and OLS reported as sensitivity.
%   * Bland-Altman sensitivity with time point 2 included (BA_ALT = 2:9): LOA and Deming slope (lambda re-estimated).
%   * B0 for the PRF conversion comes from the DICOM ImagingFrequency (gre_params.m), not nominal 3 T.
%
% REVISION 2026-11-23 (supporting statistics, macros appended only)
%   * S5b same-subset LOOCV without time point 2 (GRE vs 2DVSMSI on the same time points), near/remote pooled.
%   * S5c far-to-near calibration with the real GRE PRF signal at ROIs I/IV (literature alpha, no probe) -> ROIs II/III.
%   * S5d ROI II position sensitivity (3x3 window shifted towards the rod end).
%   * S5e single-measurement PFA precision from the back-to-back 2b blocks (data/rev_2b_blocks.mat,
%     code/supporting_analyses/rev_2b_precision.m) and block 1 as a second PFA calibration.
%
% Sections
%   S0  configuration and data loading
%   S1  per-ROI calibration and PRF                    (tempfitT2old / PRFreconnew)
%   S1b per-ROI slopes, CV, SOS common slope, cluster bootstrap CI, paired PFA-Full difference
%   S2  LOOCV + per-ROI error metrics                  -> Table 1 / Table S1
%   S3  Figure 4  calibration curves (PRF, PFA-SOS, Full-SOS)
%       Figure S1 MIP calibration curves
%   S4  Figure 5  dT maps (NaN outside valid mask)     + Figure S2 absolute T
%   S5  Figure S3 Bland-Altman on dT (bias + LOA only)
%   S5b-e same-subset LOOCV, far-to-near calibration, ROI II position, 2b precision
%   S6  text blocks and .mat of the results
%   S7  LaTeX macros (numbers.tex)
%
% Requires on the path: tempfitT2old.m, PRFreconnew.m
% Created 2026-08-13; rewritten 2026-09-03 for the native-grid pipeline.

clear; close all;
project = project_paths();
cd(project.code);

%% ========================================================================
%  S0  CONFIGURATION AND DATA LOADING
%% ========================================================================
OUT_DIR = project.results;
FIG_DIR = project.figures;   % \graphicspath of the TBME manuscript

NB          = 20000;         % bootstrap replicates
SEED        = 20260830;      % bootstrap seed (fixed)
EDGE_W      = 10;            % border width for background-noise estimation
NOISE_FACTOR= 3;             % signal below 3x background is masked
TP_CMP      = 3;             % time point of the dT maps (Fig. 5, S2, S4), display thresholds and SNR
BA_TPS      = 3:9;           % time points pooled in the Bland-Altman comparison (Fig. S3) and
                             % averaged as the dT reference of the maps and the Bland-Altman.
                             % Time point 1 is not used: its GRE preceded its 2DVSMSI scans by
                             % ~5 min.  Time point 2 is not used: its GRE phase scan (series
                             % 0017) deviates strongly from the cooling trend of the other scans.
REF_ALT     = 3:9;           % alternative dT2/ds reference used only by the invariance check (S1a)
N128        = 128;  N100 = 100;

CLIM_T  = [35 50];           % absolute temperature maps (supplementary), degC
CLIM_DT = [-5 10];           % dT maps (main Figure 5), degC; set to match the
                             % range of the absolute-temperature colour bar
                             % (T - Tref spans roughly -3 to +8 degC here).

% Purple box (SNR / Bland-Altman region), defined on the 128 grid as in the
% original analysis: rows 51:90, cols 45:50.
BOX_ROWS_128 = 51:90;  BOX_COLS_128 = 45:50;

% Native-grid rod and ROI coordinates are resolved by analysis_geometry.

% ---- data -----------------------------------------------------------------
S1f = load(fullfile(project.data,'T2MAP_and_img4i_native100.mat'), 'T2MAP','T2MAPsos','img4i','img4isos','reconstruction_frame','fit_seconds');
R   = load(fullfile(project.data,'ROI_and_tempr_resolved.mat'));
S3f = load(fullfile(project.data,'YX.mat'), 'X','Y','yangaib');
S4f = load(fullfile(project.data,'phaseMAPr_and_magmap.mat'), 'phaseMAPr','magmap');
ND  = load(fullfile(project.data,'noise_dependence_summary.mat'));
if ~isfield(ND,'tp') || ND.tp ~= TP_CMP
    warning('stats:SnrTimePoint',['noise_dependence_summary.mat was not computed at time point %d. SNR entries are set to NaN. ' ...
        'Fix: run rev_pseudoreplica.m (PR_TP=%d), rev_pseudoreplica_analysis.m and rev_write_dependence_summary.m, then rerun this script.'],TP_CMP,TP_CMP);
    ND.snr_replica = nan(1,4);
end

T2MAP    = double(S1f.T2MAP);       % MIP bin combination, 100x100x9x3
T2MAPsos = double(S1f.T2MAPsos);    % SOS bin combination
img4isos = double(S1f.img4isos);    % first-echo SOS magnitude
img4i = double(S1f.img4i);         % first fitted-echo MIP magnitude
phaseMAPr= S4f.phaseMAPr;           % 128x128x9, GRE phase
magmap   = double(S4f.magmap);      % 128x128, GRE magnitude

temp1 = R.temp1;  temp3 = R.temp3;  temprp1 = R.temprp1;     % 9x4, degC; mean fluoroptic-log temperature over each scan's data-acquisition window, 51 s clock offset corrected (build_probe_temperatures.m)
oil1 = R.oil1; oil2 = R.oil2; oil3 = R.oil3; oil4 = R.oil4;   % 128 grid, PRF drift ref
yangaib128 = logical(S3f.yangaib);                             % 128 grid, display mask

n_tp = size(T2MAPsos,3);  n_roi = 4;  n_vox = 9;

% ---- ROI coordinates in the verified reconstruction frame ---------------
geometry = analysis_geometry(R,S1f.reconstruction_frame);
C2D = geometry.pfa_roi;
C2D_full = geometry.full_roi;
ROD_R=geometry.rod_100(1:2); ROD_C=geometry.rod_100(3:4);
PX_MM = geometry.px_100;                  % native voxel [readout PE], mm
%          % 4x2 [row col], native grid
CPRF = geometry.prf_roi;                  % 4x2 [row col], 128 grid
roi_2d = geometry.pfa_series;
roi_full = geometry.full_series;
%        % {9x2} x4, for tempfitT2old
roi_prf = R.ROI_PRF_128_series;            % {9x2} x4, for PRFreconnew
assert(size(T2MAPsos,3)==size(phaseMAPr,3) && size(temp1,1)==size(T2MAPsos,3) && size(temp3,1)==size(T2MAPsos,3) ...
    && size(temprp1,1)==size(T2MAPsos,3) && size(roi_2d{1},1)==size(T2MAPsos,3) && size(roi_prf{1},1)==size(T2MAPsos,3), ...
    'stats:TimePoints','Inconsistent number of time points between maps, temperatures and ROI series. Fix: rerun reconstruct_all.m and build_probe_temperatures.m.');
fprintf('ROI 2DVSMSI (native 100): %s\n', mat2str(C2D));
fprintf('ROI Full (native 100)  : %s\n', mat2str(C2D_full));
fprintf('ROI PRF     (128 grid)  : %s\n', mat2str(CPRF));

% Descriptive magnitude SNR at time point TP_CMP, before interpolation.  The only
% GRE magnitude image is series 0016 (acquisition 2); it is used for every GRE
% magnitude quantity (no GRE magnitude exists for the other acquisitions).
% Use the same physical reference rectangle and an empty upper-left corner.
snr_rows100 = unique(round(gre_to_native(BOX_ROWS_128, 1, geometry)));
snr_cols100 = unique(round(gre_to_native(BOX_COLS_128, 2, geometry)));
fprintf('reference rectangle on the native grid: rows %d-%d, cols %d-%d\n', ...
    snr_rows100(1), snr_rows100(end), snr_cols100(1), snr_cols100(end));
apparent_snr = [magnitude_snr(magmap,BOX_ROWS_128,BOX_COLS_128,1:13), ...
    magnitude_snr(img4isos(:,:,TP_CMP,1),snr_rows100,snr_cols100,1:10), ...
    magnitude_snr(img4isos(:,:,TP_CMP,3),snr_rows100,snr_cols100,1:10), ...
    magnitude_snr(img4i(:,:,TP_CMP,1),snr_rows100,snr_cols100,1:10), ...
    magnitude_snr(img4i(:,:,TP_CMP,3),snr_rows100,snr_cols100,1:10)];   % legacy patch-SD definition
% Rayleigh-corrected GRE background (sigma = mean/sqrt(pi/2)); the four 13x13
% corners give the range because the GRE background contains structured ghosts.
gre_sig = mean(double(magmap(BOX_ROWS_128,BOX_COLS_128)),'all');
crn = {1:13, N128-12:N128};
gre_snr_corners = zeros(1,4); q = 0;
for a_ = 1:2, for b_ = 1:2
    q = q + 1; bgp = double(magmap(crn{a_},crn{b_}));
    gre_snr_corners(q) = gre_sig/(mean(bgp(:))/sqrt(pi/2));
end, end
snr_revised = [gre_snr_corners(1), ND.snr_replica];   % [GRE PFA-SOS Full-SOS PFA-MIP Full-MIP]

% ---- oil mask for display: yangaib from YX.mat covers the full bottle outline,
%      while oil1..oil4 from ROI_and_tempr_resolved.mat only mark the bright cores
%      (224 vs 5291 voxels on the 128 grid).  Use yangaib for mask exclusion so
%      no bottle edge survives into the temperature maps; keep oil1..oil4 for PRF
%      drift correction (they provide four independent reference regions).
[ri,ci]    = ndgrid(1:N100, 1:N100);
r128 = min(max(round(native_to_gre(ri, 1, geometry)),1),N128);
c128 = min(max(round(native_to_gre(ci, 2, geometry)),1),N128);
oil2d_pfa  = yangaib128(sub2ind([N128 N128], r128, c128));
oil2d_pfa = circshift(oil2d_pfa,geometry.oil_from_legacy);
oil2d_full = oil2d_pfa; % Both reconstructions now use centered-FFT coordinates.

% ---- distances to rod surface, for Table 1 --------------------------------
dist_mm = zeros(1,n_roi);
for k = 1:n_roi
    dist_mm(k) = rod_dist(C2D(k,:), ROD_R(1),ROD_R(2),ROD_C(1),ROD_C(2), PX_MM, k);
end
fprintf('distance to rod surface (mm): %s  -> paper: %s\n', ...
    mat2str(dist_mm,3), mat2str(round(dist_mm)));

%% ========================================================================
%  S1  PER-ROI CALIBRATION
%% ========================================================================
fprintf('\n---- S1  per-ROI calibration -----------------------------------\n');
[dMNS_p_sos, SDS_p_sos, ~, ~, ~, m2_p_sos] = tempfitT2old(squeeze(T2MAPsos(:,:,:,1)), roi_2d{:}, temp1);
[dMNS_f_sos, SDS_f_sos, ~, ~, ~, m2_f_sos] = tempfitT2old(squeeze(T2MAPsos(:,:,:,3)), roi_full{:}, temp3);
[dMNS_p_mip, SDS_p_mip, ~, ~, ~, m2_p_mip] = tempfitT2old(squeeze(T2MAP(:,:,:,1)),    roi_2d{:}, temp1);
[dMNS_f_mip, SDS_f_mip, ~, ~, ~, m2_f_mip] = tempfitT2old(squeeze(T2MAP(:,:,:,3)),    roi_full{:}, temp3);

[MNSp, SDSp, ~, ~, ~, alpham] = PRFreconnew(phaseMAPr, roi_prf{:}, temprp1, oil1,oil2,oil3,oil4);
ALPHA_PAPER = -alpham;      % positive, ds per degC (Eq. 7 sign convention)
MNS_DS      = -MNSp;        % ROI-mean ds deviation

prt_slopes('PFA-SOS ', m2_p_sos);  prt_slopes('Full-SOS', m2_f_sos);
prt_slopes('PFA-MIP ', m2_p_mip);  prt_slopes('Full-MIP', m2_f_mip);
fprintf('alpha PRF   [%s] x1e-8/degC   remote (I,IV) mean %.3f\n', ...
    sprintf('%.3f ',ALPHA_PAPER*1e8), mean(ALPHA_PAPER([1 4]))*1e8);

fit_pfa  = roifit(dMNS_p_sos, temp1);   fit_full = roifit(dMNS_f_sos, temp3);
fit_pmip = roifit(dMNS_p_mip, temp1);   fit_fmip = roifit(dMNS_f_mip, temp3);
fit_prf  = roifit(MNS_DS,     temprp1);

%% ========================================================================
%  S1a  REFERENCE INVARIANCE: dT2 / ds referenced to all time points vs REF_ALT
%% ========================================================================
% The per-voxel reference only adds a constant to every ROI series, which the
% ROI intercepts absorb.  Recompute the deviations with the reference taken over
% REF_ALT through the same code paths and compare slopes, homogeneity F and
% LOOCV errors.
fprintf('\n---- S1a reference invariance (all %d time points vs %d-%d) ----\n', size(T2MAPsos,3), REF_ALT(1), REF_ALT(end));
[dalt_p_sos,~,~,~,~,m2alt_p_sos] = tempfitT2old(squeeze(T2MAPsos(:,:,:,1)), roi_2d{:},   temp1, REF_ALT);
[dalt_f_sos,~,~,~,~,m2alt_f_sos] = tempfitT2old(squeeze(T2MAPsos(:,:,:,3)), roi_full{:}, temp3, REF_ALT);
[dalt_p_mip,~,~,~,~,m2alt_p_mip] = tempfitT2old(squeeze(T2MAP(:,:,:,1)),    roi_2d{:},   temp1, REF_ALT);
[dalt_f_mip,~,~,~,~,m2alt_f_mip] = tempfitT2old(squeeze(T2MAP(:,:,:,3)),    roi_full{:}, temp3, REF_ALT);
[MNSp_alt,~,~,~,~,alpham_alt]    = PRFreconnew(phaseMAPr, roi_prf{:}, temprp1, oil1,oil2,oil3,oil4, REF_ALT);
alt_cfg = {-MNSp_alt, temprp1; dalt_p_sos, temp1; dalt_f_sos, temp3; dalt_p_mip, temp1; dalt_f_mip, temp3};
ref_cfg = {MNS_DS, temprp1; dMNS_p_sos, temp1; dMNS_f_sos, temp3; dMNS_p_mip, temp1; dMNS_f_mip, temp3};
inv_names = {'GRE','PFA-SOS','Full-SOS','PFA-MIP','Full-MIP'};
inv_slope = zeros(1,5); inv_loocv = zeros(1,5); inv_offset = zeros(1,5);
for m = 1:5
    [sla,mca] = common_fit(alt_cfg{m,1}, alt_cfg{m,2});
    [slb,mcb] = common_fit(ref_cfg{m,1}, ref_cfg{m,2});
    ea = loocv_err(alt_cfg{m,1}, alt_cfg{m,2});  eb = loocv_err(ref_cfg{m,1}, ref_cfg{m,2});
    inv_slope(m) = max(abs([sla-slb, mca-mcb]));
    inv_loocv(m) = max(abs(ea-eb),[],'all');
    inv_offset(m) = max(abs(mean(alt_cfg{m,1},1)-mean(ref_cfg{m,1},1)));   % the constant that differs (should be nonzero)
    fprintf('%-9s max|slope diff| %.2e  max|LOOCV error diff| %.2e degC   (ROI-mean offset between references up to %.3g)\n', ...
        inv_names{m}, inv_slope(m), inv_loocv(m), inv_offset(m));
end
assert(max(inv_slope) < 1e-9 && max(inv_loocv) < 1e-7, 'stats:RefInvariance', ...
    'Reference choice changed slopes or LOOCV. Fix: check tempfitT2old/PRFreconnew reft handling.');
ref_check = struct('names',{inv_names},'ref_alt',REF_ALT,'max_slope_diff',inv_slope,'max_loocv_diff',inv_loocv,'ref_offset',inv_offset);

%% ========================================================================
%  S1b  PER-ROI SLOPES, CV, SOS COMMON SLOPE, CLUSTER BOOTSTRAP CI
%% ========================================================================
% Descriptive only.  No F test, no power analysis, no p-value gate: those assume
% independent ROI residuals, which this design does not support (see pair_r below).
% Documented rule (no threshold): the common coefficient is reported for the SOS
% combination (primary); MIP slopes are reported per ROI, because their between-ROI CV
% was much larger (printed below).  The bootstrap resamples whole time points, so the four
% ROIs of a time point stay together and ROI independence is not assumed.
fprintf('\n---- S1b per-ROI slopes / CV / SOS common slope / bootstrap ----\n');
cal = struct('name',{'PFA-SOS','Full-SOS','PFA-MIP','Full-MIP'}, ...
             'Y',{dMNS_p_sos, dMNS_f_sos, dMNS_p_mip, dMNS_f_mip}, ...
             'T',{temp1, temp3, temp1, temp3});
POOL = struct([]);
bootstrap_indices = [];
for c = 1:numel(cal)
    [sl,mc,Rres] = common_fit(cal(c).Y, cal(c).T);
    rf = roifit(cal(c).Y, cal(c).T);
    POOL(c).name = cal(c).name;
    POOL(c).slopes = sl;   POOL(c).cv = 100*std(sl)/mean(sl);
    POOL(c).r2 = rf.r2;    POOL(c).pooled = mc;
    POOL(c).common_reported = any(strcmp(cal(c).name,{'PFA-SOS','Full-SOS'}));   % documented descriptive rule
    if POOL(c).common_reported
        [~,ci,bootstrap_indices] = pooled_boot(cal(c).Y, cal(c).T, NB, SEED);
        POOL(c).ci = ci;  POOL(c).relhw = 100*(ci(2)-ci(1))/2/mc;
    else
        POOL(c).ci = [NaN NaN];  POOL(c).relhw = NaN;
    end
    fprintf('%-9s per-ROI slopes [%s]  R^2 [%s]  CV %.1f%%  common %.4f  CI [%.4f, %.4f]  %s\n', ...
        POOL(c).name, sprintf('%.3f ',sl), sprintf('%.3f ',rf.r2), POOL(c).cv, mc, POOL(c).ci, ...
        tern(POOL(c).common_reported,'common slope reported (SOS)','per-ROI only (MIP)'));
end
% pairwise residual correlation across ROIs (SOS only): why the bootstrap resamples whole time points
for c = 1:2
    [~,~,Rres] = common_fit(cal(c).Y, cal(c).T);
    cm = corrcoef(Rres);
    v = [cm(1,2) cm(1,3) cm(1,4) cm(2,3) cm(2,4) cm(3,4)];
    POOL(c).pair_r = v;
    fprintf('%-9s pairwise residual r: I-II %+.2f I-III %+.2f I-IV %+.2f II-III %+.2f II-IV %+.2f III-IV %+.2f\n', ...
        POOL(c).name, v);
end
RES_CORR_MAX = max([POOL(1).pair_r POOL(2).pair_r]);
% Paired cluster bootstrap of the PFA - Full common-slope difference (same resampled time points for both).
[PAIR_DIFF, PAIR_CI, PAIR_NVALID] = paired_boot_diff(cal(1).Y, cal(1).T, cal(2).Y, cal(2).T, NB, SEED);
assert(abs(PAIR_DIFF - (POOL(1).pooled-POOL(2).pooled)) < 1e-12, 'stats:PairDiff', 'Paired difference does not equal the difference of the common slopes. Fix: check paired_boot_diff.');
fprintf('PFA-SOS minus Full-SOS common slope %.4f, paired cluster-bootstrap 95%% CI [%.4f, %.4f] (%d valid draws)\n', ...
    PAIR_DIFF, PAIR_CI, PAIR_NVALID);

%% ========================================================================
%  S2  LOOCV + PER-ROI ERROR METRICS                    -> Table 1 / S1
%% ========================================================================
fprintf('\n---- S2  LOOCV (leave one time point out) ----------------------\n');
loocv_cfg = { ...
    {MNS_DS,     temprp1, 'GRE'     }, ...
    {dMNS_p_sos, temp1,   'PFA-SOS' }, ...
    {dMNS_f_sos, temp3,   'Full-SOS'}, ...
    {dMNS_p_mip, temp1,   'PFA-MIP' }, ...
    {dMNS_f_mip, temp3,   'Full-MIP'} };
n_m = numel(loocv_cfg);
loocv_errors = zeros(n_tp,n_roi,n_m);
loocv_rmse = zeros(n_m,n_roi); loocv_mae = zeros(n_m,n_roi);
loocv_bias = zeros(n_m,n_roi); loocv_max = zeros(n_m,n_roi);
for m = 1:n_m
    MNS = loocv_cfg{m}{1};  Tp = loocv_cfg{m}{2};
    for k = 1:n_roi
        e = zeros(n_tp,1);
        for lo = 1:n_tp
            tr = setdiff(1:n_tp, lo);
            xt = Tp(tr,k) - mean(Tp(tr,k));  yt = MNS(tr,k) - mean(MNS(tr,k));
            bt = sum(xt.*yt)/sum(xt.^2);
            T_pred = (MNS(lo,k) - mean(MNS(tr,k)))/bt + mean(Tp(tr,k));
            e(lo) = T_pred - Tp(lo,k);
        end
        loocv_rmse(m,k) = sqrt(mean(e.^2));
        loocv_mae(m,k)  = mean(abs(e));
        loocv_bias(m,k) = mean(e);
        loocv_max(m,k)  = max(abs(e));
        loocv_errors(:,k,m) = e;
    end
end
fprintf('%-9s %8s %8s %8s %8s   RMS\n','LOOCV RMSE','ROI I','ROI II','ROI III','ROI IV');
for m = 1:n_m
    fprintf('%-9s %8.2f %8.2f %8.2f %8.2f  %6.2f\n', loocv_cfg{m}{3}, loocv_rmse(m,:), sqrt(mean(loocv_rmse(m,:).^2)));
end

% ---- S2b  secondary metrics: common-coefficient LOOCV, bootstrap intervals, in-sample SD ----
% (i) Common-coefficient LOOCV: each fold refits the model that produced the temperature maps (one common slope plus
%     ROI intercepts) on the other time points.  Computed for the SOS variants and, to show why it is not used, MIP.
% (ii) 95% cluster-bootstrap intervals for the overall LOOCV RMSE (whole time points of the held-out errors resampled,
%     NB draws, rng(SEED), randi(n,n,1) per draw: same convention as the slope bootstrap).
% (iii) In-sample SD of (calibration-fit temperature - probe temperature) per ROI, N-1 normalisation, the quantity
%     Weber et al. 2017 call thermometry uncertainty.  SOS: the common-slope fit used for the maps; GRE: per-ROI fit.
fprintf('\n---- S2b common-coefficient LOOCV, bootstrap intervals, in-sample SD ----\n');
common_idx = [2 3 4 5];                       % PFA-SOS, Full-SOS, PFA-MIP, Full-MIP
common_errors = zeros(n_tp,n_roi,numel(common_idx));
common_rmse = zeros(numel(common_idx),n_roi);
for q = 1:numel(common_idx)
    common_errors(:,:,q) = loocv_common(loocv_cfg{common_idx(q)}{1}, loocv_cfg{common_idx(q)}{2});
    common_rmse(q,:) = sqrt(mean(common_errors(:,:,q).^2,1));
    fprintf('%-9s common-coef LOOCV RMSE %s  overall %.2f   (per-ROI-calibration overall %.2f)\n', loocv_cfg{common_idx(q)}{3}, ...
        sprintf('%.2f ',common_rmse(q,:)), sqrt(mean(common_rmse(q,:).^2)), sqrt(mean(loocv_rmse(common_idx(q),:).^2)));
end
% intervals: GRE remote (ROIs I, IV), GRE near (II, III), four variants (all ROIs), four common-coefficient variants
boot_sets = {loocv_errors(:,[1 4],1), loocv_errors(:,[2 3],1), loocv_errors(:,:,2), loocv_errors(:,:,3), loocv_errors(:,:,4), loocv_errors(:,:,5), ...
             common_errors(:,:,1), common_errors(:,:,2), common_errors(:,:,3), common_errors(:,:,4)};
boot_names = {'GreRemote','GreNear','PfaSos','FullSos','PfaMip','FullMip','CommonPfaSos','CommonFullSos','CommonPfaMip','CommonFullMip'};
boot_rmse = zeros(numel(boot_sets),1);  boot_ci = zeros(numel(boot_sets),2);
for q = 1:numel(boot_sets)
    boot_rmse(q) = sqrt(mean(boot_sets{q}.^2,'all'));
    boot_ci(q,:) = boot_rmse_ci(boot_sets{q}, NB, SEED);
    fprintf('overall LOOCV RMSE %-14s %.2f  95%% cluster-bootstrap CI [%.2f, %.2f]\n', boot_names{q}, boot_rmse(q), boot_ci(q,:));
end
% in-sample SD of fit-derived minus probe temperature
sd_fit_roi = zeros(n_m,n_roi);         % per-ROI fit (every variant)
sd_fit_common = zeros(n_m,n_roi);      % common-slope fit (all but GRE)
for m = 1:n_m
    Yq = loocv_cfg{m}{1}; Tq = loocv_cfg{m}{2};
    [slq,mcq] = common_fit(Yq,Tq);
    for k = 1:n_roi
        yc = Yq(:,k)-mean(Yq(:,k));  tc = Tq(:,k)-mean(Tq(:,k));
        sd_fit_roi(m,k)    = std(yc/slq(k) - tc);
        sd_fit_common(m,k) = std(yc/mcq - tc);
    end
    fprintf('%-9s in-sample SD (per-ROI fit) %s | (common fit) %s\n', loocv_cfg{m}{3}, sprintf('%.2f ',sd_fit_roi(m,:)), sprintf('%.2f ',sd_fit_common(m,:)));
end

%% ========================================================================
%  S3  FIGURE 4 (main: PRF / PFA-SOS / Full-SOS)  and  FIGURE S1 (MIP)
%% ========================================================================
fprintf('\n---- S3  Figure 4 / S1 ------------------------------------------\n');
% Error bars: +/-1 SD (N-1 normalisation, as returned by tempfitT2old and
% PRFreconnew) across the nine native-grid voxels of each ROI. Slope CIs use
% the separate cluster bootstrap.
fprintf('error bars: +/-1 SD across the %d ROI voxels\n', n_vox);

rows_main = { ...
  {MNS_DS,     SDSp,      ALPHA_PAPER, temprp1, '\Deltas (10^{-8})',   'GRE',                  fit_prf,  1e8}, ...
  {dMNS_p_sos, SDS_p_sos, m2_p_sos,    temp1,   '\DeltaT_{2,app} (ms)','2DVSMSI-PFA (SOS)',    fit_pfa,  1  }, ...
  {dMNS_f_sos, SDS_f_sos, m2_f_sos,    temp3,   '\DeltaT_{2,app} (ms)','2DVSMSI-Full (SOS)',   fit_full, 1  } };
fig4 = calib_figure(rows_main, n_roi, POOL(1:2));
exportgraphics(fig4, fullfile(FIG_DIR,'Fig4_calibration_SOS.png'), 'Resolution',300,'BackgroundColor','white');

rows_supp = { ...
  {dMNS_p_mip, SDS_p_mip, m2_p_mip, temp1, '\DeltaT_{2,app} (ms)','2DVSMSI-PFA (MIP)',  fit_pmip, 1}, ...
  {dMNS_f_mip, SDS_f_mip, m2_f_mip, temp3, '\DeltaT_{2,app} (ms)','2DVSMSI-Full (MIP)', fit_fmip, 1} };
figS1 = calib_figure(rows_supp, n_roi, []);
exportgraphics(figS1, fullfile(FIG_DIR,'FigS1_calibration_MIP.png'), 'Resolution',300,'BackgroundColor','white');

%% ========================================================================
%  S4  FIGURE 5  dT maps (NaN outside mask)  +  FIGURE S2 absolute T
%% ========================================================================
fprintf('\n---- S4  Figure 5 / S2 ------------------------------------------\n');
% Display threshold = NOISE_FACTOR x the maximum of the EDGE_W-pixel edge band
% of each method's own magnitude at time point TP_CMP (16 ms SOS for 2DVSMSI).
% The GRE magnitude image is series 0016 (acquisition 2): the only GRE magnitude
% that exists, used for the GRE threshold whatever TP_CMP is.
noise_pfa  = edge_noise(img4isos(:,:,TP_CMP,1), EDGE_W);
noise_full = edge_noise(img4isos(:,:,TP_CMP,3), EDGE_W);
noise_gre  = edge_noise(magmap,            EDGE_W);
ref_sig = [mean(img4isos(snr_rows100,snr_cols100,TP_CMP,1),'all'), ...
    mean(img4isos(snr_rows100,snr_cols100,TP_CMP,3),'all'), gre_sig];
thr_disp = NOISE_FACTOR*[noise_pfa noise_full noise_gre];
fprintf('display threshold %dx edge max: PFA %.4g (%.1f%% of ref. rectangle) | Full %.4g (%.1f%%) | GRE %.4g (%.1f%%)\n', ...
    NOISE_FACTOR, [thr_disp; 100*thr_disp./ref_sig]);

Tref_pfa  = mean(temp1(BA_TPS,:),'all');  Tref_full = mean(temp3(BA_TPS,:),'all');  Tref_prf = mean(temprp1(BA_TPS,:),'all');
m_pfa  = POOL(1).pooled;   m_full = POOL(2).pooled;     % pooled SOS slopes
a_prf  = mean(ALPHA_PAPER([1 4]));                       % remote-ROI alpha

[T_pfa,  mk_pfa,  dT_pfa ] = t2_maps(T2MAPsos(:,:,:,1), img4isos(:,:,TP_CMP,1), m_pfa,  Tref_pfa,  oil2d_pfa,  noise_pfa,  NOISE_FACTOR, TP_CMP, BA_TPS);
[T_full, mk_full, dT_full] = t2_maps(T2MAPsos(:,:,:,3), img4isos(:,:,TP_CMP,3), m_full, Tref_full, oil2d_full, noise_full, NOISE_FACTOR, TP_CMP, BA_TPS);
[T_prf,  mk_prf,  dT_prf ] = prf_maps(phaseMAPr, magmap, a_prf, Tref_prf, noise_gre, NOISE_FACTOR, yangaib128, oil1,oil2,oil3,oil4, TP_CMP, BA_TPS);
fprintf('valid voxels: PFA %d | Full %d (native 100) | PRF %d (128)\n', nnz(mk_pfa), nnz(mk_full), nnz(mk_prf));

% Resample the native maps onto the GRE grid by physical position (centred
% grids, anisotropic native voxel). Used for display and for the paired
% Bland-Altman comparison only; calibration uses the native grid.
% NaN is preserved (invalid voxels stay blank).
dT_pfa_d  = resample_to_gre(dT_pfa,  N128, geometry);   dT_full_d = resample_to_gre(dT_full, N128, geometry);
T_pfa_d   = resample_to_gre(T_pfa,   N128, geometry);   T_full_d  = resample_to_gre(T_full,  N128, geometry);

% POCS now returns the standard centered-FFT image origin. Both variants
% are already in the same frame; applying the old display shift would shift
% PFA a second time. Interpolation changes grid spacing, not the origin.
dT_pfa_r = dT_pfa_d;
dT_full_r = dT_full_d;
rods = {geometry.rod_128([3 4 1 2]), ...
    geometry.rod_resampled_128([3 4 1 2]), geometry.rod_resampled_128([3 4 1 2])};

% Purple box (the Bland-Altman reference region), drawn as a hollow outline
% on every panel.  On the 128 grid it is rows 51:90, cols 45:50 (the original
% 240-voxel region mentioned in the manuscript).
BOX = [BOX_COLS_128(1)-0.5, BOX_ROWS_128(1)-0.5, ...
       numel(BOX_COLS_128), numel(BOX_ROWS_128)];

% Maps are drawn on a plain black axis with no magnitude underlay.  An earlier
% revision underlaid the magnitude image, but the only thing it ever added were
% the four oil bottles and the sub-threshold rim, neither of which is part of
% the temperature result.  Voxels outside the mask are NaN and therefore read
% as background here, unlike the S4b snapshots where they take the clim floor.
fig5 = figure('Color','w','Position',[60 60 1150 400]);
tiledlayout(1,3,'TileSpacing','compact','Padding','compact');
maps  = {dT_prf, dT_pfa_r, dT_full_r};
names = {'GRE','2DVSMSI-PFA (SOS)','2DVSMSI-Full (SOS)'};
for k = 1:3
    ax = nexttile(k);
    h = imagesc(maps{k});  axis image; hold on;
    set(h,'AlphaData',~isnan(maps{k}),'AlphaDataMapping','none');
    set(ax,'Color','k','XTick',[],'YTick',[],'XColor','none','YColor','none');
    colormap(ax, parula);  clim(ax, CLIM_DT);
    rod = rods{k};
    patch([rod(1) rod(1) rod(2) rod(2)], [rod(3) rod(4) rod(4) rod(3)], ...
        'r','FaceColor','red','EdgeColor','none');
    rectangle('Position',BOX,'EdgeColor','m','LineWidth',0.8,'LineStyle','-');
    hold off; title(names{k},'FontSize',10);
end
cb = colorbar(gca,'eastoutside'); ylabel(cb,'\DeltaT (\circC)','FontSize',10);
exportgraphics(fig5, fullfile(FIG_DIR,'Fig5_dT_maps.png'), 'Resolution',300,'BackgroundColor','white');

figS2 = figure('Color','w','Position',[60 60 1150 400]);
tiledlayout(1,3,'TileSpacing','compact','Padding','compact');
% Both maps already use the corrected centered-FFT origin.
T_pfa_dr = T_pfa_d;
maps  = {T_prf, T_pfa_dr, T_full_d};
for k = 1:3
    ax = nexttile(k);
    h = imagesc(maps{k});  axis image; hold on;
    set(h,'AlphaData',~isnan(maps{k}),'AlphaDataMapping','none');
    set(ax,'Color','k','XTick',[],'YTick',[],'XColor','none','YColor','none');
    colormap(ax, parula);  clim(ax, CLIM_T);
    rod = rods{k};
    patch([rod(1) rod(1) rod(2) rod(2)], [rod(3) rod(4) rod(4) rod(3)], ...
        'r','FaceColor','red','EdgeColor','none');
    rectangle('Position',BOX,'EdgeColor','m','LineWidth',0.8,'LineStyle','-');
    hold off;
    title(names{k},'FontSize',10);
end
cb = colorbar(gca,'eastoutside'); ylabel(cb,'Temperature (\circC)','FontSize',10);
exportgraphics(figS2, fullfile(FIG_DIR,'FigS2_absolute_T_maps.png'), 'Resolution',300,'BackgroundColor','white');

%% ========================================================================
%  S4b  FIGURE 4 subpanels B/C/D — grayscale snapshot maps with colorbar
%       prfmap.png  = Ds map  (GRE,   native 128, time point TP_CMP)
%       dt2mapB.png = DT2 map (PFA,   native 100, time point TP_CMP)
%       dt2mapA.png = DT2 map (Full,  native 100, time point TP_CMP)
%
% Conventions for these three panels (author decision 2026-09-11):
%   * NATIVE grids, no display interpolation and no PE-line registration.
%     Each panel shows the map exactly as reconstructed on its own grid.
%   * Masked voxels are ASSIGNED the clim minimum (-3 ms / -3e-8) rather than
%     NaN, reproducing the manuscript convention ("we assigned a dT2 value of
%     -3 ms to voxels with signal intensities below three times the background
%     noise").  No NaN survives into these images.
%   * Masks are the ORIGINAL per-modality intensity thresholds returned by
%     t2_maps / prf_maps (mk_pfa / mk_full / mk_prf): magnitude at or above
%     3x the edge noise with the oil references removed, plus the T2 fit
%     validity terms on the T2 branch.  One rule across all three panels.
%% ========================================================================
fprintf('\n---- S4b  Figure 4 B/C/D (snapshot maps) -----------------------\n');

% Rod polygon on each panel's own grid, from the physical size (8 cm x 1 cm):
%   GRE  128 grid, 2.00 mm/voxel -> cols 62-67, rows 51-91
%   PFA and Full 100 grid, 2.56 mm/voxel -> cols 49-52, rows 40-70
rod100 = geometry.rod_100([3 4 1 2]);
rod128 = geometry.rod_128([3 4 1 2]);
ROD_S4B = {rod128, rod100, rod100};

% Colour ranges: dT2 symmetric ±3 ms; Ds ±6×10⁻⁸ (dTE = 2.46 ms; ±3×10⁻⁸ saturated).  Masked voxels sit exactly at the clim
% minimum, so the gray colormap renders them at its dark end and they read as
% the -3 reference value; the figure background itself stays white.
CLIM_DS  = [-6 6]*1e-8;   % Ds range: ±6×10⁻⁸ (was ±3×10⁻⁸ before the dTE fix 5->2.46 ms, which scales Ds by 5/2.46)
CLIM_DT2 = [-3 3];        % DT2 range: ±3 ms (matches manuscript text)
BG_DS    = CLIM_DS(1);    % masked-voxel fill value, Ds  map
BG_DT2   = CLIM_DT2(1);   % masked-voxel fill value, DT2 map

% --- Ds map (prfmap.png) ------------------------------------------------
ph_vol = double(phaseMAPr);
ph_uw  = permute(unwrap(permute(ph_vol,[3 1 2])),[2 3 1]);
dph    = ph_uw - mean(ph_uw(:,:,BA_TPS),3);
oil_all = logical(max(cat(3,oil1,oil2,oil3,oil4),[],3));
id_oil = find(oil_all);
for kk = 1:size(dph,3)
    sl = dph(:,:,kk); dph(:,:,kk) = sl - mean(sl(id_oil),'all');
end
gp_=gre_params(); lam_=gp_.lam;  % dTE = 2.46 ms, see gre_params.m
ds_map = -dph(:,:,TP_CMP)/lam_;
% Same mask as prf_maps (3x background noise, oil references excluded), and
% masked voxels take the clim minimum instead of NaN.
ds_disp = ds_map;
ds_disp(~mk_prf) = BG_DS;

figPRF = figure('Color','w','Position',[50 50 560 520]);
axP = axes(figPRF,'Units','normalized','Position',[0.08 0.06 0.68 0.90]);
imagesc(axP, ds_disp); axis(axP,'image','off');
colormap(axP, gray); clim(axP, CLIM_DS);
cbP = colorbar(axP,'eastoutside'); cbP.FontSize = 9;
cbP.Label.String = '\Delta s (\times10^{-8})';
cbP.Ticks = (-6:3:6)*1e-8;
cbP.TickLabels = {'-6','-3','0','3','6'};
hold(axP,'on');
rod_ = ROD_S4B{1};
patch(axP,[rod_(1) rod_(1) rod_(2) rod_(2)],[rod_(3) rod_(4) rod_(4) rod_(3)],...
    'r','FaceColor','red','EdgeColor','none');
hold(axP,'off');
exportgraphics(figPRF, fullfile(FIG_DIR,'prfmap.png'), 'Resolution',300,'BackgroundColor','white');
close(figPRF); fprintf('wrote prfmap.png\n');

% --- DT2 map PFA (dt2mapB.png) ------------------------------------------
% Native 100 grid, no interpolation and no PE registration: this panel shows
% the PFA reconstruction in its own frame.
dT2_pfa_tp1 = (T2MAPsos(:,:,TP_CMP,1) - mean(T2MAPsos(:,:,BA_TPS,1),3));   % ms
dT2_pfa_tp1(~mk_pfa) = BG_DT2;

figDTB = figure('Color','w','Position',[50 50 560 520]);
axB = axes(figDTB,'Units','normalized','Position',[0.08 0.06 0.68 0.90]);
imagesc(axB, dT2_pfa_tp1); axis(axB,'image','off');
colormap(axB, gray); clim(axB, CLIM_DT2);
cbB = colorbar(axB,'eastoutside'); cbB.FontSize = 9;
cbB.Label.String = '\DeltaT_2 (ms)';
cbB.Ticks = -3:1.5:3; cbB.TickLabels = {'-3','-1.5','0','1.5','3'};
hold(axB,'on');
rod_ = ROD_S4B{2};
patch(axB,[rod_(1) rod_(1) rod_(2) rod_(2)],[rod_(3) rod_(4) rod_(4) rod_(3)],...
    'r','FaceColor','red','EdgeColor','none');
hold(axB,'off');
exportgraphics(figDTB, fullfile(FIG_DIR,'dt2mapB.png'), 'Resolution',300,'BackgroundColor','white');
close(figDTB); fprintf('wrote dt2mapB.png\n');

% --- DT2 map Full (dt2mapA.png) -----------------------------------------
% Native 100 grid, no interpolation; the Full frame is used as reconstructed.
dT2_full_tp1 = (T2MAPsos(:,:,TP_CMP,3) - mean(T2MAPsos(:,:,BA_TPS,3),3));   % ms
dT2_full_tp1(~mk_full) = BG_DT2;

figDTA = figure('Color','w','Position',[50 50 560 520]);
axA = axes(figDTA,'Units','normalized','Position',[0.08 0.06 0.68 0.90]);
imagesc(axA, dT2_full_tp1); axis(axA,'image','off');
colormap(axA, gray); clim(axA, CLIM_DT2);
cbA = colorbar(axA,'eastoutside'); cbA.FontSize = 9;
cbA.Label.String = '\DeltaT_2 (ms)';
cbA.Ticks = -3:1.5:3; cbA.TickLabels = {'-3','-1.5','0','1.5','3'};
hold(axA,'on');
rod_ = ROD_S4B{3};
patch(axA,[rod_(1) rod_(1) rod_(2) rod_(2)],[rod_(3) rod_(4) rod_(4) rod_(3)],...
    'r','FaceColor','red','EdgeColor','none');
hold(axA,'off');
exportgraphics(figDTA, fullfile(FIG_DIR,'dt2mapA.png'), 'Resolution',300,'BackgroundColor','white');
close(figDTA); fprintf('wrote dt2mapA.png\n');

%% ========================================================================
%  S5  FIGURE S3  Bland-Altman on dT, native grid, bias + LOA only
%% ========================================================================
fprintf('\n---- S5  Bland-Altman (native grid, purple box) -----------------\n');
% dT referenced to each voxel's own temporal mean; the unknown spatial
% baseline cancels.  PRF is the reference, so it stays on its native 128
% grid; the 2DVSMSI dT maps are interpolated up to 128 for the comparison
% (display/paired-comparison only - no quantitative statistic is derived
% per interpolated voxel).  dT_pfa_r / dT_full_r are the PE-registered maps
% computed above, so both are on the PRF frame.
%
% Purple box on the 128 grid: rows 51-90, cols 45-50 (the original box).

[rr128,cc128] = ndgrid(1:N128, 1:N128);
box = false(N128);  box(BOX_ROWS_128, BOX_COLS_128) = true;
% masks: valid in the method AND in the reference box
% native valid -> 128 grid, nearest native voxel at the same physical position
r100 = round(gre_to_native(rr128, 1, geometry));  c100 = round(gre_to_native(cc128, 2, geometry));
inside = r100>=1 & r100<=N100 & c100>=1 & c100<=N100;
mk_pfa_128 = false(N128);  mk_full_128 = false(N128);
mk_pfa_128(inside)  = mk_pfa(sub2ind([N100 N100], r100(inside), c100(inside)));
mk_full_128(inside) = mk_full(sub2ind([N100 N100], r100(inside), c100(inside)));

sel = box & mk_pfa_128 & mk_full_128 & mk_prf;
% Pool the voxel pairs of every time point in BA_TPS; per-time-point bias and
% LOA are kept for the supplementary table.  The masks do not depend on time.
a = []; b = []; p = [];
BA_tp = zeros(numel(BA_TPS),6);      % [PFA bias lo hi, Full bias lo hi]
for it = 1:numel(BA_TPS)
    tp = BA_TPS(it);
    [~,~,da] = t2_maps(T2MAPsos(:,:,:,1), img4isos(:,:,TP_CMP,1), m_pfa,  Tref_pfa,  oil2d_pfa,  noise_pfa,  NOISE_FACTOR, tp, BA_TPS);
    [~,~,db] = t2_maps(T2MAPsos(:,:,:,3), img4isos(:,:,TP_CMP,3), m_full, Tref_full, oil2d_full, noise_full, NOISE_FACTOR, tp, BA_TPS);
    [~,~,dp] = prf_maps(phaseMAPr, magmap, a_prf, Tref_prf, noise_gre, NOISE_FACTOR, yangaib128, oil1,oil2,oil3,oil4, tp, BA_TPS);
    da = resample_to_gre(da, N128, geometry);  db = resample_to_gre(db, N128, geometry);
    at = da(sel);  bt = db(sel);  pt = dp(sel);
    s1 = bastats(at, pt);  s2 = bastats(bt, pt);
    BA_tp(it,:) = [s1.bias s1.loa_lo s1.loa_hi s2.bias s2.loa_lo s2.loa_hi];
    fprintf('  time point %d: PFA-SOS bias %+.2f LOA [%+.2f, %+.2f] | Full-SOS bias %+.2f LOA [%+.2f, %+.2f]\n', tp, BA_tp(it,:));
    a = [a; at];  b = [b; bt];  p = [p; pt]; %#ok<AGROW>
end
n_vox = nnz(sel);
n_pairs = numel(a);
fprintf('Bland-Altman on the 128 grid: purple box rows %d-%d cols %d-%d = %d voxels; %d valid x %d time points = %d pairs\n', ...
    BOX_ROWS_128(1),BOX_ROWS_128(end),BOX_COLS_128(1),BOX_COLS_128(end), nnz(box), n_vox, numel(BA_TPS), n_pairs);
BA1 = bastats(a, p);  BA2 = bastats(b, p);
fprintf('PFA-SOS  vs PRF: bias %+.2f  LOA [%+.2f, %+.2f]  RMSE %.3f\n', BA1.bias, BA1.loa_lo, BA1.loa_hi, BA1.rmse);
fprintf('Full-SOS vs PRF: bias %+.2f  LOA [%+.2f, %+.2f]  RMSE %.3f\n', BA2.bias, BA2.loa_lo, BA2.loa_hi, BA2.rmse);
% The pooled bias is zero by construction (dT referenced to the mean over
% BA_TPS), so proportional agreement is summarised by the Deming slope of
% 2DVSMSI dT on GRE dT, with a CI from resampling whole time points.
% Error-variance ratio (user decision 2026-11-17, option A): the three methods measure the same
% temperature change in the same voxels, so each method's error variance is estimated with Grubbs'
% three-instrument estimator (assumes mutually independent errors), and the Deming slope uses
% lambda = var(err 2DVSMSI) / var(err GRE), held at the full-data estimate in every bootstrap draw.
% Sensitivity: lambda = 1 (orthogonal regression) and lambda -> Inf (ordinary least squares, GRE error-free).
nt_ba = numel(BA_TPS);
[DEM_errvar, DEM_lam] = grubbs_ratio(a, b, p);     % errvar = [PFA Full GRE] (degC^2); lam = [PFA Full]
[DEM_pfa,  DEM_pfa_ci ] = deming_tp(reshape(a,n_vox,nt_ba), reshape(p,n_vox,nt_ba), NB, SEED, DEM_lam(1));
[DEM_full, DEM_full_ci] = deming_tp(reshape(b,n_vox,nt_ba), reshape(p,n_vox,nt_ba), NB, SEED, DEM_lam(2));
[DEMeq_pfa,  DEMeq_pfa_ci ] = deming_tp(reshape(a,n_vox,nt_ba), reshape(p,n_vox,nt_ba), NB, SEED, 1);
[DEMeq_full, DEMeq_full_ci] = deming_tp(reshape(b,n_vox,nt_ba), reshape(p,n_vox,nt_ba), NB, SEED, 1);
[DEMols_pfa,  DEMols_pfa_ci ] = deming_tp(reshape(a,n_vox,nt_ba), reshape(p,n_vox,nt_ba), NB, SEED, Inf);
[DEMols_full, DEMols_full_ci] = deming_tp(reshape(b,n_vox,nt_ba), reshape(p,n_vox,nt_ba), NB, SEED, Inf);
fprintf('Grubbs error variances (degC^2): PFA %.3f Full %.3f GRE %.3f; lambda PFA %.2f Full %.2f\n', DEM_errvar, DEM_lam);
fprintf('Deming slope vs PRF dT (estimated lambda): PFA-SOS %.3f [%.3f, %.3f] | Full-SOS %.3f [%.3f, %.3f] (time-point bootstrap)\n', ...
    DEM_pfa, DEM_pfa_ci, DEM_full, DEM_full_ci);
fprintf('  lambda = 1:   PFA-SOS %.3f [%.3f, %.3f] | Full-SOS %.3f [%.3f, %.3f]\n', DEMeq_pfa, DEMeq_pfa_ci, DEMeq_full, DEMeq_full_ci);
fprintf('  lambda = Inf: PFA-SOS %.3f [%.3f, %.3f] | Full-SOS %.3f [%.3f, %.3f]\n', DEMols_pfa, DEMols_pfa_ci, DEMols_full, DEMols_full_ci);

% Sensitivity: the same comparison with time point 2 included (BA_ALT = 2:9; dT reference = mean over
% 2:9; valid-voxel masks recomputed with that reference; common slopes unchanged).  Reported in the supplement only.
BA_ALT = 2:9;
[~,mkA_p] = t2_maps(T2MAPsos(:,:,:,1), img4isos(:,:,TP_CMP,1), m_pfa,  Tref_pfa,  oil2d_pfa,  noise_pfa,  NOISE_FACTOR, TP_CMP, BA_ALT);
[~,mkA_f] = t2_maps(T2MAPsos(:,:,:,3), img4isos(:,:,TP_CMP,3), m_full, Tref_full, oil2d_full, noise_full, NOISE_FACTOR, TP_CMP, BA_ALT);
mkA_p128 = false(N128);  mkA_f128 = false(N128);
mkA_p128(inside) = mkA_p(sub2ind([N100 N100], r100(inside), c100(inside)));
mkA_f128(inside) = mkA_f(sub2ind([N100 N100], r100(inside), c100(inside)));
selA = box & mkA_p128 & mkA_f128 & mk_prf;
a_alt = []; b_alt = []; p_alt = [];
for tp = BA_ALT
    [~,~,da] = t2_maps(T2MAPsos(:,:,:,1), img4isos(:,:,TP_CMP,1), m_pfa,  Tref_pfa,  oil2d_pfa,  noise_pfa,  NOISE_FACTOR, tp, BA_ALT);
    [~,~,db] = t2_maps(T2MAPsos(:,:,:,3), img4isos(:,:,TP_CMP,3), m_full, Tref_full, oil2d_full, noise_full, NOISE_FACTOR, tp, BA_ALT);
    [~,~,dp] = prf_maps(phaseMAPr, magmap, a_prf, Tref_prf, noise_gre, NOISE_FACTOR, yangaib128, oil1,oil2,oil3,oil4, tp, BA_ALT);
    da = resample_to_gre(da, N128, geometry);  db = resample_to_gre(db, N128, geometry);
    a_alt = [a_alt; da(selA)];  b_alt = [b_alt; db(selA)];  p_alt = [p_alt; dp(selA)]; %#ok<AGROW>
end
BAalt1 = bastats(a_alt, p_alt);  BAalt2 = bastats(b_alt, p_alt);
nv_alt = nnz(selA);  nt_alt = numel(BA_ALT);
[~, DEMalt_lam] = grubbs_ratio(a_alt, b_alt, p_alt);   % lambda re-estimated for this time-point set
[DEMalt_pfa,  DEMalt_pfa_ci ] = deming_tp(reshape(a_alt,nv_alt,nt_alt), reshape(p_alt,nv_alt,nt_alt), NB, SEED, DEMalt_lam(1));
[DEMalt_full, DEMalt_full_ci] = deming_tp(reshape(b_alt,nv_alt,nt_alt), reshape(p_alt,nv_alt,nt_alt), NB, SEED, DEMalt_lam(2));
fprintf('Sensitivity, time points %d-%d (%d voxels, %d pairs): LOA PFA-SOS [%+.2f, %+.2f], Full-SOS [%+.2f, %+.2f]; Deming PFA %.3f [%.3f, %.3f], Full %.3f [%.3f, %.3f]\n', ...
    BA_ALT(1), BA_ALT(end), nv_alt, numel(a_alt), BAalt1.loa_lo, BAalt1.loa_hi, BAalt2.loa_lo, BAalt2.loa_hi, ...
    DEMalt_pfa, DEMalt_pfa_ci, DEMalt_full, DEMalt_full_ci);

figS3 = figure('Color','w','Position',[60 60 1560 440]);
tlS3 = tiledlayout(1,3,'TileSpacing','compact','Padding','compact');
% common axis ranges so the two panels are directly comparable
ba_x = [floor(min([(a+p);(b+p)]/2)*10)/10-0.2,  ceil(max([(a+p);(b+p)]/2)*10)/10+0.2];
ba_y = [min([a-p; b-p; BA1.loa_lo; BA2.loa_lo])-0.25,  max([a-p; b-p; BA1.loa_hi; BA2.loa_hi])+0.25];
% x label is the mean of the two delta-T values being compared; per-panel
% titles are omitted so both panels read as one comparison.
ba_panel(nexttile(1), (a+p)/2, a-p, BA1, [0.20 0.40 0.80], '\DeltaT_{2DVSMSI-PFA} - \DeltaT_{GRE} (\circC)', '', ba_x, ba_y);
ba_panel(nexttile(2), (b+p)/2, b-p, BA2, [0.90 0.40 0.00], '\DeltaT_{2DVSMSI-Full} - \DeltaT_{GRE} (\circC)', '', ba_x, ba_y);
% per-time-point bias and LOA, PFA and Full offset horizontally
axT = nexttile(3); hold(axT,'on');
errorbar(axT, BA_TPS-0.12, BA_tp(:,1), BA_tp(:,1)-BA_tp(:,2), BA_tp(:,3)-BA_tp(:,1), 'o', ...
    'Color',[0.20 0.40 0.80],'MarkerFaceColor',[0.20 0.40 0.80],'LineWidth',1.1,'CapSize',4);
errorbar(axT, BA_TPS+0.12, BA_tp(:,4), BA_tp(:,4)-BA_tp(:,5), BA_tp(:,6)-BA_tp(:,4), 's', ...
    'Color',[0.90 0.40 0.00],'MarkerFaceColor',[0.90 0.40 0.00],'LineWidth',1.1,'CapSize',4);
plot(axT, [BA_TPS(1)-0.5 BA_TPS(end)+0.5], [0 0], 'k:', 'LineWidth',0.5);
xlim(axT, [BA_TPS(1)-0.5 BA_TPS(end)+0.5]);  xticks(axT, BA_TPS);
ylim(axT, [min(BA_tp(:,[2 5]),[],'all')-0.3, max(BA_tp(:,[3 6]),[],'all')+0.3]);
xlabel(axT, 'Time point','FontSize',10);
ylabel(axT, 'Per-time-point bias and 95% LOA (\circC)','FontSize',10);
legend(axT, {'2DVSMSI-PFA-SOS','2DVSMSI-Full-SOS'}, 'Location','northoutside','Orientation','horizontal','FontSize',8.5,'Box','off');
set(axT,'FontSize',9,'Box','on'); hold(axT,'off');
title(tlS3, sprintf('Bland-Altman of \\DeltaT within the reference box (time points %d-%d, n = %d voxel pairs)', BA_TPS(1), BA_TPS(end), n_pairs), 'FontSize',11);
exportgraphics(figS3, fullfile(FIG_DIR,'FigS3_bland_altman.png'), 'Resolution',300,'BackgroundColor','white');

%% ========================================================================
%  S5b  SAME-SUBSET LOOCV: GRE and 2DVSMSI on the same time points
%% ========================================================================
% The GRE near-rod LOOCV error is dominated by time point 2, whose GRE phase image is the known outlier.
% For a like-for-like comparison every method is also evaluated without time point 2 (time points 1, 3-9),
% with the per-ROI LOOCV of S2 (own slope and intercept, one time point left out), and the held-out errors
% are pooled over the near-rod ROIs (II, III) and the remote ROIs (I, IV).
fprintf('\n---- S5b same-subset LOOCV (all nine vs without time point 2) ----\n');
SUB_TPS = [1 3:n_tp];
SUB = same_subset_loocv(loocv_cfg(1:3), loocv_errors(:,:,1:3), SUB_TPS);
for q_ = 1:3
    fprintf('%-9s all nine: near %.2f remote %.2f | without time point 2: per ROI [%s] near %.2f remote %.2f degC\n', SUB(q_).name, ...
        SUB(q_).near_all, SUB(q_).remote_all, sprintf('%.2f ',SUB(q_).rmse_roi), SUB(q_).near, SUB(q_).remote);
end

%% ========================================================================
%  S5c  FAR-TO-NEAR CALIBRATION WITH THE REAL GRE PRF (no probe in the calibration)
%% ========================================================================
% Calibration route that needs no probe: the GRE PRF temperature change at the remote ROIs I and IV
% (X = ds/alpha with the literature alpha = 1e-8 /degC, no probe-fitted coefficient) is the temperature
% reference.  X is linearly interpolated in time from the GRE window midpoints to each 2DVSMSI window midpoint
% (outside the GRE range: linear extrapolation from the two nearest GRE points); only the GRE scans the scheme
% allows are used.  The common slope of the ROI-mean apparent-T2 deviation on X (ROI-specific intercepts) over
% the calibration time points converts the near-rod ROIs II and III; their probes are used only for the anchor
% temperature and for scoring.  Schemes: fwd (primary; calibrate 1-5, test 6-9, anchor 1), online_now
% (calibrate 1..j for test point j), online_past (1..j-1), rev (calibrate 5-9, test 1-4, anchor 9), and
% fwd_nointerp (fwd with X taken at the GRE times, no interpolation).
fprintf('\n---- S5c far-to-near calibration with the GRE PRF signal (literature alpha) ----\n');
ALPHA_LIT = 1e-8;                                       % literature PRF coefficient, ds per degC
FAR_CAL_LAST = 5;                                       % fwd: calibrate on time points 1-5, test 6-9
t_mid = (R.scan_start_s + R.scan_end_s)/2;              % window midpoints, 9 x [GRE PFA Full], scanner clock (s)
FAR_SCHEMES = {'fwd','online_now','online_past','rev','fwd_nointerp'};
far_in = {'PFA-SOS', dMNS_p_sos, temp1, t_mid(:,2); 'Full-SOS', dMNS_f_sos, temp3, t_mid(:,3)};
FAR = cell(1,0);
for v_ = 1:2
    for s_ = 1:numel(FAR_SCHEMES)
        f_ = far_to_near(far_in{v_,1}, far_in{v_,2}, far_in{v_,3}, far_in{v_,4}, MNS_DS, t_mid(:,1), ALPHA_LIT, FAR_SCHEMES{s_}, FAR_CAL_LAST);
        FAR{end+1} = f_; %#ok<AGROW>
        fprintf('%-9s %-12s test %s: RMSE II %.2f III %.2f both %.2f | bias %+.2f %+.2f | max %.2f | slope %s\n', f_.variant, f_.scheme, ...
            mat2str(f_.test_tp), f_.rmse, f_.rmse_both, f_.bias, f_.maxe, sprintf('%.3f ',f_.slopes));
    end
end
FAR = [FAR{:}];
% GRE itself near the rod with the literature alpha, same anchor (time point 1) and test points, no calibration
FAR_GRE.test_tp = FAR_CAL_LAST+1:n_tp;
FAR_GRE.errors = (MNS_DS(FAR_GRE.test_tp,[2 3]) - MNS_DS(1,[2 3]))/ALPHA_LIT - (temprp1(FAR_GRE.test_tp,[2 3]) - temprp1(1,[2 3]));
FAR_GRE.rmse = sqrt(mean(FAR_GRE.errors.^2,1));
fprintf('GRE near the rod, literature alpha, anchor time point 1, test %s: RMSE II %.2f III %.2f\n', mat2str(FAR_GRE.test_tp), FAR_GRE.rmse);
% range of the per-ROI RMSE over the sensitivity schemes
FAR_SENS = zeros(2,2);                                  % [min max] x [PFA; Full]
for v_ = 1:2
    sel_ = strcmp({FAR.variant}, far_in{v_,1}) & ismember({FAR.scheme}, {'online_now','online_past','rev'});
    rr_ = vertcat(FAR(sel_).rmse);
    FAR_SENS(v_,:) = [min(rr_(:)) max(rr_(:))];
    fprintf('%-9s per-ROI RMSE over online_now/online_past/rev: %.2f-%.2f degC\n', far_in{v_,1}, FAR_SENS(v_,:));
end

%% ========================================================================
%  S5d  ROI II POSITION SENSITIVITY
%% ========================================================================
% ROI II is the ROI nearest to the rod end.  Its 3x3 window is shifted by drow in {0,-1,-2} (towards the rod
% end) and dcol in {-1,0,+1}; for each position: per-ROI slope against the ROI II probe, per-ROI LOOCV RMSE,
% distance to the rod surface (rod_dist, ROI II rule) and the mean first-fitted-echo (16 ms) SOS magnitude
% relative to the unshifted ROI.
fprintf('\n---- S5d ROI II position sensitivity ----\n');
POS_DR = [0 -1 -2];  POS_DC = [-1 0 1];
POS = cell(1,2);
pos_in = {'PFA-SOS', 1, C2D(2,:), temp1(:,2); 'Full-SOS', 3, C2D_full(2,:), temp3(:,2)};
for v_ = 1:2
    POS{v_} = roi_position(pos_in{v_,1}, T2MAPsos(:,:,:,pos_in{v_,2}), img4isos(:,:,:,pos_in{v_,2}), pos_in{v_,3}, pos_in{v_,4}, ...
        POS_DR, POS_DC, ROD_R, ROD_C, PX_MM);
end
POS = [POS{:}];
assert(abs(POS(1).slope(1,2)-m2_p_sos(2)) < 1e-10 && abs(POS(2).slope(1,2)-m2_f_sos(2)) < 1e-10 ...
    && abs(POS(1).rmse(1,2)-loocv_rmse(2,2)) < 1e-10 && abs(POS(2).rmse(1,2)-loocv_rmse(3,2)) < 1e-10, 'stats:RoiPosition', ...
    'The unshifted ROI II window does not reproduce the Table 1 slope/LOOCV. Fix: check roi_position against tempfitT2old.');
for v_ = 1:2
    for i_ = 1:numel(POS_DR)
        for j_ = 1:numel(POS_DC)
            fprintf('%-9s shift (%+d,%+d) centre %s: dist %.1f mm  slope %.3f  LOOCV %.2f  magnitude %.2f x  invalid %d/81\n', POS(v_).variant, ...
                POS_DR(i_), POS_DC(j_), mat2str(POS(v_).centre(:,:,i_,j_)), POS(v_).dist(i_,j_), POS(v_).slope(i_,j_), POS(v_).rmse(i_,j_), ...
                POS(v_).mag_rel(i_,j_), POS(v_).n_invalid(i_,j_));
        end
    end
end
POS_RANGE = pos_ranges(POS, POS_DR);
fprintf('shifted towards the rod (drow < 0): PFA slope %.3f-%.3f LOOCV %.2f-%.2f | Full slope %.3f-%.3f LOOCV %.2f-%.2f | distance %.1f-%.1f mm (drow -1), %.1f-%.1f mm (drow -2) | magnitude %.2f-%.2f\n', ...
    POS_RANGE.slope(1,:), POS_RANGE.rmse(1,:), POS_RANGE.slope(2,:), POS_RANGE.rmse(2,:), POS_RANGE.dist_one, POS_RANGE.dist_two, POS_RANGE.mag);

%% ========================================================================
%  S5e  SINGLE-MEASUREMENT PFA PRECISION FROM THE BACK-TO-BACK "2b" BLOCKS
%% ========================================================================
% data/rev_2b_blocks.mat (code/supporting_analyses/rev_2b_precision.m): ROI-mean apparent T2 and probe
% temperature of the two independent back-to-back PFA blocks of the nine 2025-04-25 "2b" scans.
% r = (A1 - A2)/m - (T1 - T2) with m = the PFA-SOS common slope; single-measurement SD = SD(r)/sqrt(2)
% (difference of two independent measurements), pooled over 9 scans x 4 ROIs.  Block 1 alone is also
% analysed as a second PFA calibration series with the S1b/S2 methods.
fprintf('\n---- S5e single-measurement PFA precision (2b back-to-back blocks) ----\n');
f2b_ = fullfile(project.data,'rev_2b_blocks.mat');
if isfile(f2b_)
    B2B = load(f2b_);
else
    warning('stats:Missing2b',['%s not found. The 2b precision entries are set to NaN. ' ...
        'Fix: run code/supporting_analyses/rev_2b_precision.m (needs the raw 2b files, see its header), then rerun this script.'], f2b_);
    B2B = [];
end
PREC = twob_precision(B2B, POOL(1).pooled, NB, SEED);
if PREC.available
    fprintf('single-measurement SD %.3f degC (%d pairs; per ROI %s), mean r %+.3f, SD of the difference %.3f, probe change between blocks %.2f-%.2f degC (slope %.4f)\n', ...
        PREC.sd_single, PREC.n, sprintf('%.3f ',PREC.sd_roi), PREC.mean_r, PREC.sd_diff, PREC.dT_range, PREC.m);
    fprintf('block 1 as a PFA calibration: slopes [%s] CV %.1f%% common %.4f [%.4f, %.4f]; LOOCV RMSE [%s] overall %.2f\n', ...
        sprintf('%.3f ',PREC.b1_slopes), PREC.b1_cv, PREC.b1_common, PREC.b1_ci, sprintf('%.2f ',PREC.b1_rmse), PREC.b1_rmse_overall);
end

%% ========================================================================
%  S6  TEXT BLOCKS + .mat OF THE RESULTS
%% ========================================================================
out = fullfile(OUT_DIR,'statistics_results_final.txt');
fid = fopen(out,'w');
fprintf(fid,'statistics_analysis_final.m  results  (native 100x100 pipeline, ROI_and_tempr_resolved.mat)\n');
fprintf(fid,'seed %d, %d bootstrap replicates\n', SEED, NB);
fprintf(fid,'=====================================================================\n\n');
fprintf(fid,'ROI 2DVSMSI (native 100): %s\n', mat2str(C2D));
fprintf(fid,'ROI Full (native 100)  : %s\n', mat2str(C2D_full));
fprintf(fid,'ROI PRF (128)           : %s\n', mat2str(CPRF));
fprintf(fid,'native voxel (mm)       : %.4f readout x %.2f phase encoding\n', PX_MM);
fprintf(fid,'distance to rod (mm)    : %s  (unrounded %s)\n\n', mat2str(round(dist_mm)), mat2str(dist_mm,3));
fprintf(fid,'Legacy patch-SD apparent SNR [GRE PFA-SOS Full-SOS PFA-MIP Full-MIP]: %s\n',mat2str(apparent_snr,5));
fprintf(fid,'Revised SNR [GRE PFA-SOS Full-SOS PFA-MIP Full-MIP]: %s  (2DVSMSI: pseudo-replica noise SD, %d/%d replicas; GRE: Rayleigh-corrected upper-left corner, corners %s)\n', ...
    mat2str(snr_revised,4), ND.n_replicas, mat2str(gre_snr_corners,3));

fprintf(fid,'TABLE 1  per-ROI slope (ms/degC or 1e-8/degC) and LOOCV RMSE (degC)\n');
fprintf(fid,'%-10s %8s %8s %8s %8s\n','','ROI I','ROI II','ROI III','ROI IV');
fprintf(fid,'%-10s %8.3f %8.3f %8.3f %8.3f  CV %.1f%%\n','PFA-SOS m',  m2_p_sos, POOL(1).cv);
fprintf(fid,'%-10s %8.3f %8.3f %8.3f %8.3f  CV %.1f%%\n','Full-SOS m', m2_f_sos, POOL(2).cv);
fprintf(fid,'%-10s %8.3f %8.3f %8.3f %8.3f  CV %.1f%%\n','PFA-MIP m',  m2_p_mip, POOL(3).cv);
fprintf(fid,'%-10s %8.3f %8.3f %8.3f %8.3f  CV %.1f%%\n','Full-MIP m', m2_f_mip, POOL(4).cv);
fprintf(fid,'%-10s %8.3f %8.3f %8.3f %8.3f\n','PRF alpha', ALPHA_PAPER*1e8);
for m = 1:n_m
    fprintf(fid,'%-10s %8.2f %8.2f %8.2f %8.2f  RMS %.2f\n', [loocv_cfg{m}{3} ' L'], loocv_rmse(m,:), sqrt(mean(loocv_rmse(m,:).^2)));
end
fprintf(fid,'\nSLOPES (descriptive; common slope reported for SOS only, MIP per ROI)\n');
for c = 1:numel(POOL)
    fprintf(fid,'%-9s slopes [%s] R2 [%s] CV %.1f%%  common %.4f  CI [%.4f, %.4f]  relHW %.2f%%\n', ...
        POOL(c).name, sprintf('%.3f ',POOL(c).slopes), sprintf('%.3f ',POOL(c).r2), POOL(c).cv, POOL(c).pooled, POOL(c).ci, POOL(c).relhw);
end
fprintf(fid,'max pairwise ROI residual correlation (SOS common-slope fits): %+.2f\n', RES_CORR_MAX);
fprintf(fid,'PFA-SOS minus Full-SOS common slope %.4f ms/degC, paired cluster-bootstrap 95%% CI [%.4f, %.4f] (%d valid draws)\n', PAIR_DIFF, PAIR_CI, PAIR_NVALID);
fprintf(fid,'\nCALIBRATION INTERVAL SENSITIVITY (SOS; not a combined error bound)\n');
for c = 1:2
    rel = POOL(c).relhw/100;
    for delta = [1 2 5 8]
        fprintf(fid,'%-9s dT %d  slope-CI half-width scale %.2f degC\n', POOL(c).name, delta, rel*delta);
    end
end
fprintf(fid,'\nTABLE S1  per-ROI error metrics (degC): RMSE / MAE / bias / max|e|\n');
for m = 1:n_m
    fprintf(fid,'%-9s RMSE [%s] MAE [%s] bias [%s] max [%s]\n', loocv_cfg{m}{3}, ...
        sprintf('%.2f ',loocv_rmse(m,:)), sprintf('%.2f ',loocv_mae(m,:)), ...
        sprintf('%+.2f ',loocv_bias(m,:)), sprintf('%.2f ',loocv_max(m,:)));
end
fprintf(fid,'\nBLAND-ALTMAN on dT pooled over time points %d-%d, 128 grid, purple box rows %d-%d cols %d-%d, %d voxels x %d time points = %d pairs\n', ...
    BA_TPS(1), BA_TPS(end), BOX_ROWS_128(1),BOX_ROWS_128(end),BOX_COLS_128(1),BOX_COLS_128(end), n_vox, numel(BA_TPS), n_pairs);
fprintf(fid,'PFA-SOS  vs PRF  bias %+.2f  LOA [%+.2f %+.2f]  RMSE %.3f\n', BA1.bias, BA1.loa_lo, BA1.loa_hi, BA1.rmse);
fprintf(fid,'Full-SOS vs PRF  bias %+.2f  LOA [%+.2f %+.2f]  RMSE %.3f\n', BA2.bias, BA2.loa_lo, BA2.loa_hi, BA2.rmse);
for it = 1:numel(BA_TPS)
    fprintf(fid,'  time point %d: PFA-SOS bias %+.2f LOA [%+.2f %+.2f] | Full-SOS bias %+.2f LOA [%+.2f %+.2f]\n', BA_TPS(it), BA_tp(it,:));
end
fprintf(fid,'dT reference = mean over time points %d-%d, so the pooled bias is zero by construction\n', BA_TPS(1), BA_TPS(end));
fprintf(fid,'Grubbs error variances (degC^2): PFA %.3f Full %.3f GRE %.3f; Deming lambda = var(err 2DVSMSI)/var(err GRE): PFA %.2f Full %.2f\n', DEM_errvar, DEM_lam);
fprintf(fid,'Deming slope (2DVSMSI dT on PRF dT, estimated lambda), time-point bootstrap 95%% CI: PFA-SOS %.3f [%.3f %.3f]  Full-SOS %.3f [%.3f %.3f]\n', ...
    DEM_pfa, DEM_pfa_ci, DEM_full, DEM_full_ci);
fprintf(fid,'  sensitivity lambda = 1 (orthogonal): PFA-SOS %.3f [%.3f %.3f]  Full-SOS %.3f [%.3f %.3f]\n', DEMeq_pfa, DEMeq_pfa_ci, DEMeq_full, DEMeq_full_ci);
fprintf(fid,'  sensitivity lambda = Inf (OLS, GRE error-free): PFA-SOS %.3f [%.3f %.3f]  Full-SOS %.3f [%.3f %.3f]\n', DEMols_pfa, DEMols_pfa_ci, DEMols_full, DEMols_full_ci);
lab_ = {'PFA-SOS ','Full-SOS','PFA-MIP ','Full-MIP'};
fprintf(fid,'\nerror bars drawn: +/-1 SD across the nine ROI voxels; SD range over %d time points x 4 ROIs\n', n_tp);
fprintf(fid,'  GRE      %.3f-%.3f x1e-8\n', min(SDSp(:))*1e8, max(SDSp(:))*1e8);
sd_ = {SDS_p_sos, SDS_f_sos, SDS_p_mip, SDS_f_mip};
for r_ = 1:4
    fprintf(fid,'  %s %.3f-%.3f ms\n', lab_{r_}, min(sd_{r_}(:)), max(sd_{r_}(:)));
end
fprintf(fid,'\nBLAND-ALTMAN SENSITIVITY with time point 2 included (BA_ALT = %d-%d, %d pairs; dT reference = mean over these)\n', BA_ALT(1), BA_ALT(end), numel(a_alt));
fprintf(fid,'Deming lambda re-estimated for these time points: PFA %.2f Full %.2f\n', DEMalt_lam);
fprintf(fid,'PFA-SOS  vs PRF  LOA [%+.2f %+.2f]  Deming %.3f [%.3f %.3f]\n', BAalt1.loa_lo, BAalt1.loa_hi, DEMalt_pfa, DEMalt_pfa_ci);
fprintf(fid,'Full-SOS vs PRF  LOA [%+.2f %+.2f]  Deming %.3f [%.3f %.3f]\n', BAalt2.loa_lo, BAalt2.loa_hi, DEMalt_full, DEMalt_full_ci);
fprintf(fid,'T_ref: PFA %.3f Full %.3f PRF %.3f degC\n', Tref_pfa, Tref_full, Tref_prf);
fprintf(fid,'display threshold (%d x max of %d-pixel edge band, own magnitude at time point %d; GRE magnitude = series 0016, acquisition 2, the only GRE magnitude): PFA %.4g (%.1f%% of ref. rectangle) Full %.4g (%.1f%%) GRE %.4g (%.1f%%)\n', ...
    NOISE_FACTOR, EDGE_W, TP_CMP, [thr_disp; 100*thr_disp./ref_sig]);
fprintf(fid,'valid voxels: PFA %d Full %d PRF %d\n', nnz(mk_pfa), nnz(mk_full), nnz(mk_prf));
fprintf(fid,'\nREFERENCE INVARIANCE: dT2/ds referenced to all %d time points vs time points %d-%d (max abs differences)\n', n_tp, REF_ALT(1), REF_ALT(end));
for m_ = 1:5
    fprintf(fid,'  %-9s slopes %.2e  LOOCV errors %.2e degC\n', inv_names{m_}, inv_slope(m_), inv_loocv(m_));
end
fprintf(fid,'time points used: calibration/LOOCV 1-%d; temperature-change maps, Bland-Altman and Deming %d-%d; illustrations at time point %d\n', n_tp, BA_TPS(1), BA_TPS(end), TP_CMP);
fprintf(fid,'\nSAME-SUBSET LOOCV (per-ROI calibration as Table 1; RMSE in degC pooled over the near-rod ROIs II+III and the remote ROIs I+IV)\n');
for q_ = 1:3
    fprintf(fid,'%-9s all nine: near %.2f remote %.2f | without time point 2 (%s): per ROI [%s] near %.2f remote %.2f\n', SUB(q_).name, ...
        SUB(q_).near_all, SUB(q_).remote_all, mat2str(SUB_TPS), sprintf('%.2f ',SUB(q_).rmse_roi), SUB(q_).near, SUB(q_).remote);
end
fprintf(fid,'\nFAR-TO-NEAR CALIBRATION WITH THE GRE PRF SIGNAL (no probe in the calibration; literature alpha %.1f x1e-8/degC; remote ROIs I+IV -> near ROIs II, III)\n', ALPHA_LIT*1e8);
fprintf(fid,'X = ds/alpha at ROIs I and IV, linearly interpolated (extrapolated outside the known range) from GRE window midpoints to each 2DVSMSI window midpoint; common slope with ROI intercepts over the calibration time points; error = (y_j - y_anchor)/m - (T_j - T_anchor)\n');
for i_ = 1:numel(FAR)
    fprintf(fid,'%-9s %-12s test %s anchor %d: RMSE II %.2f III %.2f both %.2f | bias %+.2f %+.2f | max %.2f | slope %s\n', FAR(i_).variant, FAR(i_).scheme, ...
        mat2str(FAR(i_).test_tp), FAR(i_).anchor, FAR(i_).rmse, FAR(i_).rmse_both, FAR(i_).bias, FAR(i_).maxe, sprintf('%.3f ',FAR(i_).slopes));
end
fprintf(fid,'GRE near the rod, literature alpha, no calibration, anchor time point 1, test %s: RMSE II %.2f III %.2f\n', mat2str(FAR_GRE.test_tp), FAR_GRE.rmse);
fprintf(fid,'per-ROI RMSE range over online_now/online_past/rev: PFA-SOS %.2f-%.2f, Full-SOS %.2f-%.2f\n', FAR_SENS(1,:), FAR_SENS(2,:));
fprintf(fid,'\nROI II POSITION SENSITIVITY (3x3 window shifted by drow towards the rod end and dcol; magnitude = mean 16 ms SOS magnitude relative to the unshifted window)\n');
for v_ = 1:2
    for i_ = 1:numel(POS_DR)
        for j_ = 1:numel(POS_DC)
            fprintf(fid,'%-9s shift (%+d,%+d) centre %s: dist %.1f mm  slope %.3f  LOOCV %.2f  magnitude %.2f  invalid %d/81\n', POS(v_).variant, ...
                POS_DR(i_), POS_DC(j_), mat2str(POS(v_).centre(:,:,i_,j_)), POS(v_).dist(i_,j_), POS(v_).slope(i_,j_), POS(v_).rmse(i_,j_), ...
                POS(v_).mag_rel(i_,j_), POS(v_).n_invalid(i_,j_));
        end
    end
end
fprintf(fid,'drow < 0: PFA-SOS slope %.3f-%.3f LOOCV %.2f-%.2f | Full-SOS slope %.3f-%.3f LOOCV %.2f-%.2f | distance %.1f-%.1f mm (drow -1), %.1f-%.1f mm (drow -2) | magnitude %.2f-%.2f\n', ...
    POS_RANGE.slope(1,:), POS_RANGE.rmse(1,:), POS_RANGE.slope(2,:), POS_RANGE.rmse(2,:), POS_RANGE.dist_one, POS_RANGE.dist_two, POS_RANGE.mag);
fprintf(fid,'\nSINGLE-MEASUREMENT PFA PRECISION (2025-04-25 2b scans, two back-to-back PFA blocks; data/rev_2b_blocks.mat)\n');
if PREC.available
    fprintf(fid,'MIDs %s; r = (A1-A2)/m - (T1-T2), m = PFA-SOS common slope %.4f\n', mat2str(B2B.mids), PREC.m);
    fprintf(fid,'single-measurement SD %.3f degC (%d pairs), per ROI [%s]; mean r %+.3f; SD of the difference %.3f; probe change between blocks %.2f-%.2f degC\n', ...
        PREC.sd_single, PREC.n, sprintf('%.3f ',PREC.sd_roi), PREC.mean_r, PREC.sd_diff, PREC.dT_range);
    fprintf(fid,'block 1 as a PFA calibration: slopes [%s] R2 [%s] CV %.1f%% common %.4f CI [%.4f, %.4f]; LOOCV RMSE [%s] overall %.2f\n', ...
        sprintf('%.3f ',PREC.b1_slopes), sprintf('%.3f ',PREC.b1_r2), PREC.b1_cv, PREC.b1_common, PREC.b1_ci, sprintf('%.2f ',PREC.b1_rmse), PREC.b1_rmse_overall);
    for t_ = 1:size(B2B.A1,1)
        fprintf(fid,'  MID %d: T block 1 [%s] block 2 [%s] degC; T2app block 1 [%s] block 2 [%s] ms; r [%s]\n', B2B.mids(t_), sprintf('%.3f ',B2B.T1(t_,:)), ...
            sprintf('%.3f ',B2B.T2(t_,:)), sprintf('%.3f ',B2B.A1(t_,:)), sprintf('%.3f ',B2B.A2(t_,:)), sprintf('%+.3f ',PREC.r(t_,:)));
    end
else
    fprintf(fid,'not available: data/rev_2b_blocks.mat missing (run code/supporting_analyses/rev_2b_precision.m)\n');
end
fclose(fid);

%% ========================================================================
%  S7  LATEX MACROS  ->  output/results/numbers.tex
%% ========================================================================
% One macro per reported number, formatted as it appears in the paper.  Names are
% letters only: \n<Quantity><Method><Roi|TimePoint>.  See the header written below.
M = struct('name',{},'val',{});
seqs = {'Gre','Pfa','Full'};  key5 = {'Gre','PfaSos','FullSos','PfaMip','FullMip'};
key4 = {'PfaSos','FullSos','PfaMip','FullMip'};  roiN = {'I','II','III','IV'};
allT = [temprp1(:); temp1(:); temp3(:)];
cool_min = (max(R.scan_end_s(:)) - min(R.scan_start_s(:)))/60;    % first scan start -> last scan end
temp_step = mean(-diff([temprp1 temp1 temp3],1,1),'all');          % mean step between consecutive time points
gp_b = gre_params();   % B0 from the DICOM ImagingFrequency (gre_params.m)
M = addm(M,'GreWinSec','%.1f',gp_b.window_s);  M = addm(M,'GreScanSec','%.0f',gp_b.scan_s);   % GRE data window (DICOM SliceMeasurementDuration) / protocol total scan time
M = addm(M,'BZero','%.2f',gp_b.B0);  M = addm(M,'FZeroMHz','%.2f',gp_b.f0_MHz);  M = addm(M,'GammaBar','%.3f',gp_b.gamma_bar_MHz_T);
M = addm(M,'TempMin','%.1f',min(allT));   M = addm(M,'TempMax','%.1f',max(allT));
M = addm(M,'TempSpan','%.1f',max(allT)-min(allT));
M = addm(M,'CoolMin','%.0f',cool_min);    M = addm(M,'TempStep','%.2f',temp_step);
M = addm(M,'NumTp','%d',n_tp);            M = adds(M,'NumTpWord',numword(n_tp));
M = adds(M,'NumTrainWord',numword(n_tp-1));
M = addm(M,'NumPredRoi','%d',n_tp*n_roi); M = addm(M,'NumBaTp','%d',numel(BA_TPS)); M = adds(M,'NumBaTpWord',numword(numel(BA_TPS)));
M = addm(M,'BaFirstTp','%d',BA_TPS(1));   M = addm(M,'BaLastTp','%d',BA_TPS(end));
M = addm(M,'NumPairs','%d',n_pairs);      M = addm(M,'NumBoxVox','%d',n_vox);
M = addm(M,'TpCmp','%d',TP_CMP);          M = adds(M,'TpCmpWord',numword(TP_CMP));
gap_s = R.gre_to_pfa_gap_s;
M = addm(M,'GapOneSec','%.0f',gap_s(1));  M = addm(M,'GapOneMin','%.1f',gap_s(1)/60);
M = addm(M,'GapOtherMinSec','%.0f',min(gap_s(2:end))); M = addm(M,'GapOtherMaxSec','%.0f',max(gap_s(2:end)));
M = addm(M,'WinRangeMax','%.2f',R.win_range_max);      M = addm(M,'WinDevMax','%.2f',R.win_dev_max);
for k = 1:4, M = addm(M,['Dist' roiN{k}],'%.0f',dist_mm(k)); end
if ~isfield(S1f,'fit_seconds')
    error('vsmsi:MissingFitSeconds', ['T2MAP_and_img4i_native100.mat has no fit_seconds, so the FitMedian/FitMin/FitMax ' ...
        'macros cannot be written. Regenerate the file with reconstruct_all.m (Part A).']);
end
fs_ = sum(S1f.fit_seconds,3);      % per acquisition file: MIP + SOS fitting seconds
M = addm(M,'FitMedian','%.1f',median(fs_(:))); M = addm(M,'FitMin','%.1f',min(fs_(:))); M = addm(M,'FitMax','%.1f',max(fs_(:)));
% Table I / S1 / S2: coefficients and per-ROI LOOCV metrics
slopes_ = {m2_p_sos, m2_f_sos, m2_p_mip, m2_f_mip};
for k = 1:4
    M = addm(M,['AlphaGre' roiN{k}],'%.3f',ALPHA_PAPER(k)*1e8);
    for c = 1:4, M = addm(M,['Slope' key4{c} roiN{k}],'%.3f',slopes_{c}(k)); end
    for m = 1:n_m
        M = addm(M,['Rmse' key5{m} roiN{k}],'%.2f',loocv_rmse(m,k));
        M = addm(M,['Mae'  key5{m} roiN{k}],'%.2f',loocv_mae(m,k));
        M = addm(M,['Bias' key5{m} roiN{k}],'%+.2f',loocv_bias(m,k));
        M = addm(M,['Maxe' key5{m} roiN{k}],'%.2f',loocv_max(m,k));
    end
end
M = addm(M,'AlphaGreRemote','%.3f',mean(ALPHA_PAPER([1 4]))*1e8);
M = addm(M,'AlphaGreRemoteTwo','%.2f',mean(ALPHA_PAPER([1 4]))*1e8);
for m = 1:n_m, M = addm(M,['RmsOverall' key5{m}],'%.2f',sqrt(mean(loocv_rmse(m,:).^2))); end
for q = 1:4
    for k = 1:4, M = addm(M,['RmseCommon' key4{q} roiN{k}],'%.2f',common_rmse(q,k)); end
    M = addm(M,['RmsOverallCommon' key4{q}],'%.2f',sqrt(mean(common_rmse(q,:).^2)));
end
for q = 1:numel(boot_names)
    M = addm(M,['RmsBoot' boot_names{q}],'%.2f',boot_rmse(q));
    M = addm(M,['RmsBootLo' boot_names{q}],'%.2f',boot_ci(q,1));  M = addm(M,['RmsBootHi' boot_names{q}],'%.2f',boot_ci(q,2));
end
% Weber-style thermometry uncertainty: in-sample SD of (fit-derived - probe temperature)
sd_fit_SOSPFA_ = sd_fit_roi(2,:);  sd_fit_SOSFULL_ = sd_fit_roi(3,:);   % per-ROI calibration fit, as in Weber et al. 2017
M = addm(M,'SdFitPfaSosMin','%.2f',min(sd_fit_SOSPFA_));  M = addm(M,'SdFitPfaSosMax','%.2f',max(sd_fit_SOSPFA_));
M = addm(M,'SdFitFullSosMin','%.2f',min(sd_fit_SOSFULL_)); M = addm(M,'SdFitFullSosMax','%.2f',max(sd_fit_SOSFULL_));
M = addm(M,'SdFitSosMin','%.2f',min([sd_fit_SOSPFA_ sd_fit_SOSFULL_]));  M = addm(M,'SdFitSosMax','%.2f',max([sd_fit_SOSPFA_ sd_fit_SOSFULL_]));
M = addm(M,'SdFitGreRemoteMin','%.2f',min(sd_fit_roi(1,[1 4])));  M = addm(M,'SdFitGreRemoteMax','%.2f',max(sd_fit_roi(1,[1 4])));
M = addm(M,'CoolRate','%.2f',(max(allT)-min(allT))/cool_min);              % degC per min, temperature span / cooling duration
M = addm(M,'CoolFiveMin','%.1f',5*(max(allT)-min(allT))/cool_min);         % change during a 5 min T1 map
sdT_ = [std(temp1,0,1) std(temp3,0,1)];                                       % per-ROI sample SD of probe temperature (2DVSMSI scans)
M = addm(M,'TempSdRoiMin','%.1f',min(sdT_));  M = addm(M,'TempSdRoiMax','%.1f',max(sdT_));
M = addm(M,'RmseMinSos','%.2f',min(loocv_rmse(2:3,:),[],'all'));  M = addm(M,'RmseMaxSos','%.2f',max(loocv_rmse(2:3,:),[],'all'));
M = addm(M,'RmseMinMip','%.2f',min(loocv_rmse(4:5,:),[],'all'));  M = addm(M,'RmseMaxMip','%.2f',max(loocv_rmse(4:5,:),[],'all'));
M = addm(M,'RmseGreRemoteMin','%.2f',min(loocv_rmse(1,[1 4])));   M = addm(M,'RmseGreRemoteMax','%.2f',max(loocv_rmse(1,[1 4])));
M = addm(M,'RmseGreNearMin','%.2f',min(loocv_rmse(1,[2 3])));     M = addm(M,'RmseGreNearMax','%.2f',max(loocv_rmse(1,[2 3])));
% between-ROI slope CV, per-ROI R^2, common slopes (SOS), paired PFA-Full difference, interval scales
for c = 1:4
    M = addm(M,['Cv' key4{c}],'%.1f',POOL(c).cv);
end
fits_ = {fit_prf, fit_pfa, fit_full, fit_pmip, fit_fmip};
for m = 1:5
    for k = 1:4, M = addm(M,['Rsq' key5{m} roiN{k}],'%.3f',fits_{m}.r2(k)); end
end
M = addm(M,'ResCorrMax','%.2f',RES_CORR_MAX);
M = addm(M,'PairDiff','%.3f',PAIR_DIFF);
M = addm(M,'PairDiffLo','%.3f',PAIR_CI(1));  M = addm(M,'PairDiffHi','%.3f',PAIR_CI(2));
for c = 1:2
    rel_ = POOL(c).relhw/100;
    M = addm(M,['Pool' key4{c}],'%.4f',POOL(c).pooled);   M = addm(M,['PoolTwo' key4{c}],'%.2f',POOL(c).pooled);
    M = addm(M,['PoolLo' key4{c}],'%.4f',POOL(c).ci(1));  M = addm(M,['PoolHi' key4{c}],'%.4f',POOL(c).ci(2));
    M = addm(M,['RelHw' key4{c}],'%.2f',POOL(c).relhw);
    M = addm(M,['UncOne' key4{c}],'%.2f',rel_*1);  M = addm(M,['UncEight' key4{c}],'%.2f',rel_*8);
end
% Bland-Altman, Deming, probe temperatures at the comparison time point
bseq = {'Pfa','Full'};  BAs = {BA1, BA2};
for s_ = 1:2
    M = addm(M,['BaLoa' bseq{s_}],'%.2f',(BAs{s_}.loa_hi-BAs{s_}.loa_lo)/2);
    M = addm(M,['BaBiasPool' bseq{s_}],'%.2f',round(BAs{s_}.bias,2)+0);
    M = addm(M,['BaLoPool' bseq{s_}],'%+.2f',BAs{s_}.loa_lo);  M = addm(M,['BaHiPool' bseq{s_}],'%+.2f',BAs{s_}.loa_hi);
    M = addm(M,['BaRms' bseq{s_}],'%.3f',BAs{s_}.rmse);
    c0 = (s_-1)*3;
    for it = 1:numel(BA_TPS)
        w_ = numword(BA_TPS(it));  w_(1) = upper(w_(1));
        M = addm(M,['BaBias' bseq{s_} 'Tp' w_],'%+.2f',BA_tp(it,c0+1));
        M = addm(M,['BaLo' bseq{s_} 'Tp' w_],'%+.2f',BA_tp(it,c0+2));
        M = addm(M,['BaHi' bseq{s_} 'Tp' w_],'%+.2f',BA_tp(it,c0+3));
    end
end
% sensitivity: time point 2 included in the Bland-Altman comparison
M = addm(M,'BaAltFirstTp','%d',BA_ALT(1));  M = addm(M,'NumBaAltPairs','%d',numel(a_alt));
M = addm(M,'BaAltLoaPfa','%.2f',(BAalt1.loa_hi-BAalt1.loa_lo)/2);  M = addm(M,'BaAltLoaFull','%.2f',(BAalt2.loa_hi-BAalt2.loa_lo)/2);
M = addm(M,'BaAltDemPfa','%.2f',DEMalt_pfa);   M = addm(M,'BaAltDemLoPfa','%.2f',DEMalt_pfa_ci(1));   M = addm(M,'BaAltDemHiPfa','%.2f',DEMalt_pfa_ci(2));
M = addm(M,'BaAltDemFull','%.2f',DEMalt_full); M = addm(M,'BaAltDemLoFull','%.2f',DEMalt_full_ci(1)); M = addm(M,'BaAltDemHiFull','%.2f',DEMalt_full_ci(2));
M = addm(M,'BaBiasMin','%+.2f',min(BA_tp(:,[1 4]),[],'all'));  M = addm(M,'BaBiasMax','%+.2f',max(BA_tp(:,[1 4]),[],'all'));
dem_ = {DEM_pfa DEM_pfa_ci; DEM_full DEM_full_ci};
demeq_ = {DEMeq_pfa DEMeq_pfa_ci; DEMeq_full DEMeq_full_ci};  demols_ = {DEMols_pfa DEMols_pfa_ci; DEMols_full DEMols_full_ci};
for s_ = 1:2
    M = addm(M,['Dem' bseq{s_}],'%.2f',dem_{s_,1});
    M = addm(M,['DemLo' bseq{s_}],'%.2f',dem_{s_,2}(1));  M = addm(M,['DemHi' bseq{s_}],'%.2f',dem_{s_,2}(2));
    M = addm(M,['DemLam' bseq{s_}],'%.1f',DEM_lam(s_));  M = addm(M,['DemErrVar' bseq{s_}],'%.2f',DEM_errvar(s_));
    M = addm(M,['DemEq' bseq{s_}],'%.2f',demeq_{s_,1});
    M = addm(M,['DemEqLo' bseq{s_}],'%.2f',demeq_{s_,2}(1));  M = addm(M,['DemEqHi' bseq{s_}],'%.2f',demeq_{s_,2}(2));
    M = addm(M,['DemOls' bseq{s_}],'%.2f',demols_{s_,1});
    M = addm(M,['DemOlsLo' bseq{s_}],'%.2f',demols_{s_,2}(1));  M = addm(M,['DemOlsHi' bseq{s_}],'%.2f',demols_{s_,2}(2));
end
M = addm(M,'DemErrVarGre','%.2f',DEM_errvar(3));
Tall_ = {temprp1, temp1, temp3};  Tref_ = [Tref_prf Tref_full Tref_full];  Tref_(2) = Tref_pfa;
for s_ = 1:3
    for k = 1:4, M = addm(M,['Tcmp' seqs{s_} roiN{k}],'%.1f',Tall_{s_}(TP_CMP,k)); end
    M = addm(M,['TcmpMin' seqs{s_}],'%.1f',min(Tall_{s_}(TP_CMP,:)));  M = addm(M,['TcmpMax' seqs{s_}],'%.1f',max(Tall_{s_}(TP_CMP,:)));
    M = addm(M,['Tref' seqs{s_}],'%.1f',Tref_(s_));
end
thr_pct = 100*thr_disp./ref_sig;                 % [PFA Full GRE]
M = addm(M,'ThrPctPfa','%.0f',thr_pct(1));  M = addm(M,'ThrPctFull','%.0f',thr_pct(2));  M = addm(M,'ThrPctGre','%.0f',thr_pct(3));
% ---- added 2026-11-23: same-subset LOOCV, far-to-near calibration, ROI II position, 2b precision
% Same-subset LOOCV.  GRE per-ROI values without time point 2 are GreSensRmse* (numbers_gre_sensitivity.tex) and the GRE
% all-nine pooled values are RmsBootGreNear / RmsBootGreRemote, so they are not repeated here.
for q_ = 2:3
    for k = 1:4, M = addm(M,['RmseNoTwo' key5{q_} roiN{k}],'%.2f',SUB(q_).rmse_roi(k)); end
    M = addm(M,['RmsNear' key5{q_}],'%.2f',SUB(q_).near_all);  M = addm(M,['RmsRemote' key5{q_}],'%.2f',SUB(q_).remote_all);
end
for q_ = 1:3
    M = addm(M,['RmsNearNoTwo' key5{q_}],'%.2f',SUB(q_).near);  M = addm(M,['RmsRemoteNoTwo' key5{q_}],'%.2f',SUB(q_).remote);
end
% Far-to-near calibration with the GRE PRF signal
far_key_ = struct('fwd','Fwd','online_now','OnlineNow','online_past','OnlinePast','rev','Rev','fwd_nointerp','NoInterp');
far_var_ = struct('PFA_SOS','Pfa','Full_SOS','Full');
M = addm(M,'FarAlphaLit','%.1f',ALPHA_LIT*1e8);
M = addm(M,'FarCalLastTp','%d',FAR_CAL_LAST);  M = addm(M,'FarTestFirstTp','%d',FAR_CAL_LAST+1);
for i_ = 1:numel(FAR)
    nm_ = [far_var_.(strrep(FAR(i_).variant,'-','_')) far_key_.(FAR(i_).scheme)];
    M = addm(M,['FarRmse' nm_ 'II'],'%.2f',FAR(i_).rmse(1));  M = addm(M,['FarRmse' nm_ 'III'],'%.2f',FAR(i_).rmse(2));
    M = addm(M,['FarRmse' nm_ 'Both'],'%.2f',FAR(i_).rmse_both);
    if strcmp(FAR(i_).scheme,'fwd')
        M = addm(M,['FarBias' nm_ 'II'],'%+.2f',FAR(i_).bias(1));  M = addm(M,['FarBias' nm_ 'III'],'%+.2f',FAR(i_).bias(2));
        M = addm(M,['FarMaxe' nm_],'%.2f',FAR(i_).maxe);           M = addm(M,['FarSlope' nm_],'%.3f',FAR(i_).slopes(1));
    end
end
M = addm(M,'FarGreRmseII','%.2f',FAR_GRE.rmse(1));  M = addm(M,'FarGreRmseIII','%.2f',FAR_GRE.rmse(2));
M = addm(M,'FarSensMinPfa','%.2f',FAR_SENS(1,1));   M = addm(M,'FarSensMaxPfa','%.2f',FAR_SENS(1,2));
M = addm(M,'FarSensMinFull','%.2f',FAR_SENS(2,1));  M = addm(M,'FarSensMaxFull','%.2f',FAR_SENS(2,2));
% ROI II position sensitivity (unshifted values = SlopePfaSosII, RmsePfaSosII, ... in Table 1)
for v_ = 1:2
    M = addm(M,['PosSlopeMin' bseq{v_}],'%.3f',POS_RANGE.slope(v_,1));  M = addm(M,['PosSlopeMax' bseq{v_}],'%.3f',POS_RANGE.slope(v_,2));
    M = addm(M,['PosRmseMin' bseq{v_}],'%.2f',POS_RANGE.rmse(v_,1));    M = addm(M,['PosRmseMax' bseq{v_}],'%.2f',POS_RANGE.rmse(v_,2));
    M = addm(M,['PosDistZero' bseq{v_}],'%.1f',POS(v_).dist(1,2));
end
M = addm(M,'PosDistMinOne','%.1f',POS_RANGE.dist_one(1));  M = addm(M,'PosDistMaxOne','%.1f',POS_RANGE.dist_one(2));
M = addm(M,'PosDistMinTwo','%.1f',POS_RANGE.dist_two(1));  M = addm(M,'PosDistMaxTwo','%.1f',POS_RANGE.dist_two(2));
M = addm(M,'PosMagMin','%.2f',POS_RANGE.mag(1));  M = addm(M,'PosMagMax','%.2f',POS_RANGE.mag(2));
% 2b single-measurement precision and block-1 calibration (NaN when data/rev_2b_blocks.mat is missing)
prec_ = {'PrecSd','%.2f',PREC.sd_single; 'PrecSdThree','%.3f',PREC.sd_single; 'PrecSdDiff','%.2f',PREC.sd_diff; ...
         'PrecMeanRes','%+.3f',PREC.mean_r; 'PrecProbeChangeMin','%.2f',PREC.dT_range(1); 'PrecProbeChangeMax','%.2f',PREC.dT_range(2); ...
         'PrecNumPairs','%d',PREC.n; 'BlockOnePool','%.4f',PREC.b1_common; 'BlockOnePoolTwo','%.2f',PREC.b1_common; ...
         'BlockOnePoolLo','%.4f',PREC.b1_ci(1); 'BlockOnePoolHi','%.4f',PREC.b1_ci(2); 'BlockOneCv','%.1f',PREC.b1_cv; ...
         'BlockOneRmsOverall','%.2f',PREC.b1_rmse_overall};
for k = 1:4
    prec_(end+1,:) = {['PrecSd' roiN{k}],'%.3f',PREC.sd_roi(k)}; %#ok<SAGROW>
    prec_(end+1,:) = {['BlockOneSlope' roiN{k}],'%.3f',PREC.b1_slopes(k)}; %#ok<SAGROW>
    prec_(end+1,:) = {['BlockOneRmse' roiN{k}],'%.2f',PREC.b1_rmse(k)}; %#ok<SAGROW>
end
for i_ = 1:size(prec_,1)
    if PREC.available, M = addm(M,prec_{i_,1},prec_{i_,2},prec_{i_,3}); else, M = adds(M,prec_{i_,1},'NaN'); end
end
write_numbers(fullfile(OUT_DIR,'numbers.tex'), M);
fprintf('wrote numbers.tex (%d macros)\n', numel(M));

save(fullfile(OUT_DIR,'statistics_results_final.mat'), ...
    'C2D','C2D_full','CPRF','dist_mm','geometry','apparent_snr','snr_revised','gre_snr_corners','thr_disp', ...
    'm2_p_sos','m2_f_sos','m2_p_mip','m2_f_mip','ALPHA_PAPER', ...
    'POOL','loocv_rmse','loocv_mae','loocv_bias','loocv_max','loocv_cfg','loocv_errors','bootstrap_indices', ...
    'fit_pfa','fit_full','fit_pmip','fit_fmip','fit_prf', ...
    'BAalt1','BAalt2','DEMalt_pfa','DEMalt_pfa_ci','DEMalt_full','DEMalt_full_ci','a_alt','b_alt','p_alt','BA_ALT','PAIR_DIFF','PAIR_CI','PAIR_NVALID','RES_CORR_MAX', ...
    'BA1','BA2','n_pairs','BOX_ROWS_128','BOX_COLS_128', ...
    'T_pfa','T_full','T_prf','dT_pfa','dT_full','dT_prf', ...
    'mk_pfa','mk_full','mk_prf','Tref_pfa','Tref_full','Tref_prf','SEED','NB','TP_CMP','BA_TPS','BA_tp','n_vox','DEM_pfa','DEM_pfa_ci','DEM_full','DEM_full_ci', ...
    'a','b','p','DEM_errvar','DEM_lam','DEMeq_pfa','DEMeq_pfa_ci','DEMeq_full','DEMeq_full_ci','DEMols_pfa','DEMols_pfa_ci','DEMols_full','DEMols_full_ci','DEMalt_lam', ...
    'common_rmse','common_errors','boot_names','boot_rmse','boot_ci','sd_fit_roi','sd_fit_common','thr_pct','ref_check','REF_ALT','n_tp','M','cool_min','temp_step','fs_', ...
    'SUB','SUB_TPS','FAR','FAR_GRE','FAR_SENS','ALPHA_LIT','FAR_CAL_LAST','POS','POS_DR','POS_DC','POS_RANGE','PREC');
fprintf('\nwrote TBME figures to\n  %s\nand statistics_results_final.txt / .mat to\n  %s\n', FIG_DIR, OUT_DIR);

%% ========================================================================
%  helpers
%% ========================================================================
function s = tern(c,a,b), if c, s=a; else, s=b; end, end

function M = addm(M,name,fmt,val)
% Append macro \n<name> with the value formatted by fmt.
assert(~isempty(regexp(name,'^[A-Za-z]+$','once')), 'stats:MacroName', 'Macro name "%s" must be letters only. Fix: spell digits out.', name);
assert(~any(strcmp({M.name},name)), 'stats:MacroDup', 'Duplicate macro name %s. Fix: make the name unique.', name);
assert(isfinite(val), 'stats:MacroNaN', 'Macro %s is not finite. Fix: check the quantity it reports.', name);
M(end+1) = struct('name',name,'val',sprintf(fmt,val));
end

function M = adds(M,name,str)
assert(~any(strcmp({M.name},name)), 'stats:MacroDup', 'Duplicate macro name %s.', name);
M(end+1) = struct('name',name,'val',str);
end

function w = numword(n)
words = {'one','two','three','four','five','six','seven','eight','nine','ten','eleven','twelve'};
assert(n>=1 && n<=numel(words), 'stats:NumWord', 'numword supports 1-12.');
w = words{n};
end

function write_numbers(file, M)
fid = fopen(file,'w');
for tries = 1:5, if fid > 0, break; end, pause(1); fid = fopen(file,'w'); end   % cloud-sync file lock
assert(fid > 0, 'stats:Fopen', 'Cannot open %s for writing. Fix: close the file or wait for the sync client.', file);
fprintf(fid,'%% numbers.tex -- GENERATED by code/statistics_analysis_final.m; do not edit by hand.\n');
fprintf(fid,'%% Every number reported in the TBME main text and supplement is one macro here.\n');
fprintf(fid,'%% Naming: \\n<Quantity><Method><Roi|Tp...>, letters only.\n');
fprintf(fid,'%%   Method : Gre | PfaSos | FullSos | PfaMip | FullMip ; Seq (Bland-Altman, temperatures): Pfa | Full | Gre\n');
fprintf(fid,'%%   Roi    : I | II | III | IV ;  Tp<Word>: time point spelled out (TpThree)\n');
fprintf(fid,'%%   Quantity: AlphaGre (1e-8/degC), Slope (ms/degC), Rmse/Mae/Bias/Maxe (degC, Bias signed), RmsOverall,\n');
fprintf(fid,'%%     Cv (between-ROI slope CV), TempSdRoiMin/Max (per-ROI sample SD of probe temperature), Rsq (per-ROI R^2), ResCorrMax, PairDiff/PairDiffLo/PairDiffHi (PFA-SOS minus Full-SOS common slope, paired cluster bootstrap), Pool/PoolLo/PoolHi/PoolTwo (common slope, CI), RelHw,\n');
fprintf(fid,'%%     UncOne/UncEight (CI half-width scale at dT = 1 / 8 degC), BaLoa/BaBiasPool/BaLoPool/BaHiPool/BaRms (pooled),\n');
fprintf(fid,'%%     BaBias/BaLo/BaHi<Seq>Tp<Word> (Table S3), BaBiasMin/Max, Dem/DemLo/DemHi (Deming, Grubbs-estimated lambda), DemLam/DemErrVar<Seq>/DemErrVarGre, DemEq*/DemOls* (lambda = 1 / OLS sensitivity), Tcmp<Seq><Roi> (probe temperature at the\n');
fprintf(fid,'%%     comparison time point), TcmpMin/Max<Seq>, Tref, ThrPct<Seq>, Fit*, Temp*, Num*, Gap*, Win*, Dist<Roi>.\n');
fprintf(fid,'%%   Signs: Bias*/Ba*Lo/Ba*Hi/BaBiasMin/Max carry an explicit sign as plain ASCII text; wrap in $...$ for a minus sign.\n');
fprintf(fid,'%%   Use \\nName{} (not \\nName followed by a letter) in running text.\n');
fprintf(fid,'%% SNR, GRE residual SDs and GRE sensitivity are written by the supporting scripts (numbers_*.tex) and included at the end.\n');
for i = 1:numel(M)
    fprintf(fid,'\\newcommand{\\n%s}{%s}\n', M(i).name, M(i).val);
end
fprintf(fid,'\\InputIfFileExists{results/numbers_snr.tex}{}{}\n');
fprintf(fid,'\\InputIfFileExists{results/numbers_gre_residual.tex}{}{}\n');
fprintf(fid,'\\InputIfFileExists{results/numbers_gre_sensitivity.tex}{}{}\n');
fclose(fid);
end

function e = loocv_err(MNS, Tp)
% Leave-one-time-point-out errors (estimated - probe), n_tp x n_roi.
n_tp_ = size(Tp,1);  e = zeros(size(Tp));
for k = 1:size(Tp,2)
    for lo = 1:n_tp_
        tr = setdiff(1:n_tp_, lo);
        xt = Tp(tr,k) - mean(Tp(tr,k));  yt = MNS(tr,k) - mean(MNS(tr,k));
        bt = sum(xt.*yt)/sum(xt.^2);
        e(lo,k) = (MNS(lo,k) - mean(MNS(tr,k)))/bt + mean(Tp(tr,k)) - Tp(lo,k);
    end
end
end

function e = loocv_common(Y, X)
% Leave-one-time-point-out errors (estimated - probe) with the common-slope model (one slope for all ROIs,
% ROI-specific intercepts) refitted on the training time points of every fold.  n_tp x n_roi.
n_ = size(X,1);  e = zeros(size(X));
for lo = 1:n_
    tr = setdiff(1:n_, lo);
    xt = X(tr,:) - mean(X(tr,:),1);  yt = Y(tr,:) - mean(Y(tr,:),1);
    bc = sum(xt.*yt,'all')/sum(xt.^2,'all');
    e(lo,:) = (Y(lo,:) - mean(Y(tr,:),1))/bc + mean(X(tr,:),1) - X(lo,:);
end
end

function ci = boot_rmse_ci(E, NB, seed)
% 95% percentile interval of the overall RMSE of held-out errors E (n_tp x n_roi), resampling whole time points.
rng(seed); n = size(E,1); bs = zeros(NB,1);
for b = 1:NB
    ix = randi(n,n,1);
    bs(b) = sqrt(mean(E(ix,:).^2,'all'));
end
ci = quantile(bs,[0.025 0.975]);
end

function value = magnitude_snr(im,rows,cols,background_indices)
signal = im(rows,cols); background = im(background_indices,background_indices);
sigma = std(background(:));
if sigma <= 0, value = NaN; else, value = mean(signal(:))/sigma; end
end

function prt_slopes(lab, m)
fprintf('m2 %-9s [%s] mean %.4f  CV %.1f%%\n', lab, sprintf('%.4f ',m), mean(m), 100*std(m)/mean(m));
end

function f = roifit(Y, T)
n = size(Y,2);  f.slope = zeros(1,n);  f.r2 = zeros(1,n);
for k = 1:n
    x = T(:,k)-mean(T(:,k));  y = Y(:,k)-mean(Y(:,k));
    b = sum(x.*y)/sum(x.^2);  r = y - b*x;
    f.slope(k) = b;  f.r2(k) = 1 - sum(r.^2)/max(sum(y.^2),eps);
end
end

function [sl,mc,R] = common_fit(Y,X)
% Per-ROI slopes (sl), common slope with ROI-specific intercepts (mc) and the common-slope residuals (R, n x K).
n=size(Y,1); K=size(Y,2); sl=zeros(1,K); num=0; den=0;
for k=1:K
    xb=mean(X(:,k)); yb=mean(Y(:,k));
    sl(k)=sum((X(:,k)-xb).*(Y(:,k)-yb))/sum((X(:,k)-xb).^2);
    num=num+sum((X(:,k)-xb).*(Y(:,k)-yb)); den=den+sum((X(:,k)-xb).^2);
end
mc=num/den; R=zeros(n,K);
for k=1:K
    R(:,k)=(Y(:,k)-mean(Y(:,k)))-mc*(X(:,k)-mean(X(:,k)));
end
end

function [mc,ci,indices] = pooled_boot(Y,X,NB,seed)
% Cluster bootstrap: resample the time points as intact units carrying all
% four ROIs, so whatever cross-ROI dependence exists is preserved.
mc=anc(Y,X); rng(seed); bs=nan(NB,1); n=size(Y,1);
indices = zeros(n,NB,'uint8');
for b=1:NB
    ix=randi(n,n,1);
    indices(:,b) = uint8(ix);
    if any(std(X(ix,:))<1e-6), continue; end
    bs(b)=anc(Y(ix,:),X(ix,:));
end
bs=bs(~isnan(bs)); ci=quantile(bs,[0.025 0.975]);
end

function mc = anc(Y,X)
num=0; den=0;
for k=1:size(Y,2)
    xb=mean(X(:,k)); yb=mean(Y(:,k));
    num=num+sum((X(:,k)-xb).*(Y(:,k)-yb)); den=den+sum((X(:,k)-xb).^2);
end
mc=num/den;
end

function [d0,ci,nvalid] = paired_boot_diff(Y1,X1,Y2,X2,NB,seed)
% Paired cluster bootstrap of the difference of two common slopes (variant 1 minus variant 2).
% The SAME resampled time points (rows) are used for both variants, so the whole time point is the
% resampling unit and the pairing is preserved.  Same seed and draw rule as pooled_boot.  Reported as
% difference and percentile CI only (no p value).
d0=anc(Y1,X1)-anc(Y2,X2); rng(seed); bs=nan(NB,1); n=size(Y1,1);
for b=1:NB
    ix=randi(n,n,1);
    if any(std(X1(ix,:))<1e-6) || any(std(X2(ix,:))<1e-6), continue; end
    bs(b)=anc(Y1(ix,:),X1(ix,:))-anc(Y2(ix,:),X2(ix,:));
end
bs=bs(~isnan(bs)); nvalid=numel(bs); ci=quantile(bs,[0.025 0.975]);
end

function d = rod_dist(c,R1,R2,C1,C2,PX,k)
% Distance from ROI centre to the nearest rod surface.  ROI I/IV lie beside
% the rod (perpendicular distance); II/III lie near its ends.  The rod box
% indices are treated as surface coordinates.  PX = [row col] voxel size, mm.
r=c(1); col=c(2);
if r>=R1 && r<=R2
    dperp=min(abs(col-C1),abs(col-C2))*PX(2); if col>C1&&col<C2, dperp=0; end
else, dperp=NaN; end
corners=[R1 C1;R1 C2;R2 C1;R2 C2];
dcorner=min(hypot((corners(:,1)-r)*PX(1),(corners(:,2)-col)*PX(2)));
if col>=C1 && col<=C2
    dend=min(abs(r-R1),abs(r-R2))*PX(1); if r>R1&&r<R2, dend=0; end
else, dend=NaN; end
if k==1||k==4
    if ~isnan(dperp), d=dperp; else, d=dcorner; end
else
    cd=[dcorner dend]; cd=cd(~isnan(cd)); d=min(cd);
end
end

function fig = calib_figure(rows, n_roi, poolinfo)
nr = numel(rows);
fig = figure('Color','w','Position',[40 40 1180 260*nr+40]);
tiledlayout(nr,4,'TileSpacing','compact','Padding','compact');
for r = 1:nr
    c = rows{r};
    Y = c{1}; SD = c{2}; b = c{3}; Tp = c{4}; ylab = c{5}; nm = c{6}; ft = c{7}; sc = c{8};
    for k = 1:n_roi
        nexttile((r-1)*4 + k);
        x = Tp(:,k);  y = Y(:,k)*sc;  e = SD(:,k)*sc;
        is_gre = (sc == 1e8);
        xf = linspace(36, 45, 20);
        if is_gre
            % Time point 2: its GRE phase deviated from the cooling trend (excluded from the
            % cross-modality comparison).  Open marker; the red dashed line still fits all nine
            % time points (Table I), the grey line omits time point 2 (= GreSens* macros).
            o = true(size(x)); o(2) = false;
            errorbar(x(o), y(o), e(o), 'r.', 'MarkerSize', 11, 'LineWidth', 0.9); hold on;
            errorbar(x(2), y(2), e(2), 'o', 'Color', 'r', 'MarkerSize', 4.5, 'MarkerFaceColor', 'w', 'LineWidth', 0.9);
            p2 = polyfit(x(o), y(o), 1);
            plot(xf, polyval(p2, xf), '-', 'Color', [0.45 0.45 0.45], 'LineWidth', 1.0);
            fprintf('Fig4 GRE ROI %s: alpha without time point 2 = %.3f x1e-8/degC\n', roman(k), p2(1));
        else
            errorbar(x, y, e, 'r.', 'MarkerSize', 11, 'LineWidth', 0.9); hold on;
        end
        plot(xf, b(k)*sc*(xf-mean(x)) + mean(y), 'r--', 'LineWidth', 1.1);
        xlim([36 46]);
        if sc == 1e8 && (k == 2 || k == 3), ylim([-40 40]); elseif sc == 1e8, ylim([-6 6]); else, ylim([-3 3]); end  % Ds rows scaled x5/2.46 vs old ([-20 20]/[-3 3]) after the dTE fix
        if sc == 1e8
            txt = sprintf('\\alpha = %.2f \\times10^{-8}/\\circC\nR^2 = %.3f', b(k)*sc, ft.r2(k));
        else
            txt = sprintf('m_{2,app} = %.2f ms/\\circC\nR^2 = %.3f', b(k), ft.r2(k));
        end
        yl = ylim; text(36.3, yl(1)+0.84*diff(yl), txt, 'FontSize', 7.5, 'Interpreter','tex');
        if is_gre
            text(36.3, yl(1)+0.71*diff(yl), sprintf('without time point 2: \\alpha = %.2f', p2(1)), ...
                'FontSize', 7.5, 'Interpreter','tex', 'Color', [0.45 0.45 0.45]);
        end
        if k == 1, ylabel(ylab, 'FontSize', 10); end
        if r == nr, xlabel('Probe temperature (\circC)', 'FontSize', 10); end
        title(sprintf('%s  ROI %s', nm, roman(k)), 'FontSize', 9);
        set(gca,'FontSize',8); box on; hold off;
    end
    if ~isempty(poolinfo) && r >= 2
        pi_ = poolinfo(r-1);
        % Draw the pooled slope inside the ROI IV axis so it never collides
        % with the x-axis label or the next row of tiles.
        axk = nexttile((r-1)*4 + n_roi);
        axes(axk);
        yl = ylim;
        text(36.3, yl(1)+0.06*diff(yl), ...
            sprintf('common m_{2,app} = %.3f [%.3f, %.3f] ms/\\circC', pi_.pooled, pi_.ci), ...
            'FontSize',7.5,'VerticalAlignment','bottom','Interpreter','tex');
    end
end
end

function s = roman(k), r = {'I','II','III','IV'}; s = r{k}; end

function n = edge_noise(img, w)
img = double(img);
e = [reshape(img(1:w,:),[],1); reshape(img(end-w+1:end,:),[],1); ...
     reshape(img(:,1:w),[],1); reshape(img(:,end-w+1:end),[],1)];
n = max(e,[],'all');
end

function [T, mask, D] = t2_maps(T2vol, magvol, m2, Tref, oilmask, noise, nf, tp, reft)
% Absolute temperature (Eq.6) and dT maps; invalid voxels -> NaN.
% dT at time point tp is referenced to the mean over time points reft.
T2vol = double(T2vol);  T2ref = mean(T2vol(:,:,reft),3);
ok_all = all(T2vol > 0, 3);
ok_ref = (T2ref > 0) & (T2ref <= 200);
ok_mag = double(magvol(:,:,1)) >= nf*noise;   % magvol: 100x100 magnitude at the comparison time point
mask   = ok_all & ok_ref & ok_mag & ~oilmask;
D = (T2vol(:,:,tp) - T2ref)/m2;        % dT at time point tp, degC
T = D + Tref;
D(~mask) = NaN;  T(~mask) = NaN;
end

function [T, mask, D] = prf_maps(phvol, magmap, alpha_paper, Tref, noise, nf, oil_display, o1,o2,o3,o4, tp, reft)
gp = gre_params();  lamda = gp.lam;  % B0*gamma*dTE with dTE = 2.46 ms (gre_params.m)
ph  = permute(unwrap(permute(double(phvol),[3 1 2])),[2 3 1]);
dph = ph - mean(ph(:,:,reft),3);
id  = find(double(max(cat(3,o1,o2,o3,o4),[],3)));
for k = 1:size(dph,3)
    t = dph(:,:,k);  dph(:,:,k) = t - mean(t(id),'all');
end
ds   = -dph(:,:,tp)/lamda;
% Mask = the manuscript's own criterion: magnitude at or above nf x the
% background noise, oil bottles excluded.  oil_display (yangaib from YX.mat)
% covers the full bottle outline (5291 voxels on the 128 grid); o1..o4 from
% ROI_and_tempr_resolved.mat only mark the bright cores (224 voxels) and are
% used solely for drift correction above.  Using oil_display here removes
% bottle edge artefacts that survive when only the cores are masked.  This is
% the same intensity + oil-exclusion rule the T2 branch applies (t2_maps), so
% PRF, PFA and Full share one masking criterion; the T2 branch only adds
% fit-validity terms (T2>0, T2ref<=200) that have no PRF analogue.
mask = (double(magmap) >= nf*noise) & ~oil_display;
D = ds/alpha_paper;  T = D + Tref;
D(~mask) = NaN;  T(~mask) = NaN;
end

function out = resample_to_gre(map, N, g)
% Bilinear resampling of a native map onto the N x N GRE grid at the same
% physical positions; NaN regions and positions outside the native FOV are NaN.
[x,y]   = meshgrid(1:size(map,2), 1:size(map,1));
[xq,yq] = meshgrid(gre_to_native(1:N, 2, g), gre_to_native(1:N, 1, g));
valid = ~isnan(map);  f = map; f(~valid) = 0;
num = interp2(x,y,f,xq,yq,'linear');
den = interp2(x,y,double(valid),xq,yq,'linear');
out = num./den;  out(~(den >= 0.5)) = NaN;
end

function p = gre_to_native(i, dim, g)
% GRE index -> native index along dim (1 = readout rows, 2 = PE columns).
p = g.centre_100 + (i - g.centre_128)*g.px_128(dim)/g.px_100(dim);
end

function i = native_to_gre(p, dim, g)
% Native index -> GRE index along dim (1 = readout rows, 2 = PE columns).
i = g.centre_128 + (p - g.centre_100)*g.px_100(dim)/g.px_128(dim);
end

function S = bastats(x, y)
d = x - y;  n = numel(d);
S.bias = mean(d);  S.sd = std(d);
S.loa_lo = S.bias - 1.96*S.sd;  S.loa_hi = S.bias + 1.96*S.sd;
S.rmse = sqrt(mean(d.^2));  S.n = n;
end

function [s, ci] = deming_tp(X, P, nb, seed, lam)
% Deming slope of X on P with error-variance ratio lam = var(err X)/var(err P)
% (lam = 1: orthogonal regression; lam = Inf: ordinary least squares of X on P);
% columns are time points.  Whole time points are resampled because voxels within
% one time point are not independent; lam is held fixed in every draw.  Draws with
% fewer than three distinct time points are skipped.
if nargin < 5, lam = 1; end
s = deming(X(:), P(:), lam);
rng(seed);  nt = size(X,2);  bs = nan(nb,1);
for i = 1:nb
    c = randi(nt, 1, nt);
    if numel(unique(c)) < 3, continue; end
    bs(i) = deming(reshape(X(:,c),[],1), reshape(P(:,c),[],1), lam);
end
ci = prctile(bs(~isnan(bs)), [2.5 97.5]);
end

function s = deming(x, p, lam)
% Deming slope of x on p; lam = var(err x)/var(err p).  lam = Inf gives cov(p,x)/var(p).
if nargin < 3, lam = 1; end
sxx = var(p,1);  syy = var(x,1);  c = cov(p,x,1);  sxy = c(1,2);
if isinf(lam), s = sxy/sxx; return; end
s = (syy - lam*sxx + sqrt((syy - lam*sxx)^2 + 4*lam*sxy^2))/(2*sxy);
end

function [ev, lam] = grubbs_ratio(x, y, p)
% Grubbs' three-instrument estimator: x, y and p measure the same quantity with
% mutually independent errors, so var(err x) = var(x) - cov(x,y)*cov(x,p)/cov(y,p), etc.
% Returns ev = [var(err x) var(err y) var(err p)] and lam = ev(1:2)/ev(3).
C = cov([x(:) y(:) p(:)], 1);
ev = [C(1,1) - C(1,2)*C(1,3)/C(2,3), C(2,2) - C(1,2)*C(2,3)/C(1,3), C(3,3) - C(1,3)*C(2,3)/C(1,2)];
assert(all(ev > 0), 'grubbs_ratio: non-positive error variance; the independent-error model does not fit. Fix: report lambda = 1 and OLS only.');
lam = ev(1:2)/ev(3);
end

function ba_panel(ax, mn, df, S, col, ylab, ttl, xr, yr)
% ylab names the two quantities whose difference is on the y axis, e.g.
% '\DeltaT_{2DVSMSI-PFA} - \DeltaT_{GRE}'.  The x axis is the mean of the same
% two quantities, spelled out so the two panels are not read as sharing one
% physical quantity.
axes(ax);
fill([xr(1) xr(2) xr(2) xr(1)],[S.loa_lo S.loa_lo S.loa_hi S.loa_hi],[0.85 0.85 0.85],'FaceAlpha',0.30,'EdgeColor','none'); hold on;
plot([xr(1) xr(2)],[S.bias S.bias],'k-','LineWidth',1.4);
plot([xr(1) xr(2)],[S.loa_lo S.loa_lo],'k--','LineWidth',1.0);
plot([xr(1) xr(2)],[S.loa_hi S.loa_hi],'k--','LineWidth',1.0);
plot([xr(1) xr(2)],[0 0],'k:','LineWidth',0.5);
scatter(mn, df, 16, col, 'filled', 'MarkerFaceAlpha',0.55,'MarkerEdgeColor','none');
xlim(xr); ylim(yr);
xl = strrep(strrep(ylab,' - ',' + '),' (\circC)','');
xlabel(['(' xl ') / 2  (\circC)'],'FontSize',10);
ylabel(ylab,'FontSize',10);
if ~isempty(ttl), title(ttl,'FontSize',11,'FontWeight','bold'); end
text(xr(1)+0.08, S.bias+0.12,  sprintf('bias = %+.2f \\circC', round(S.bias,2)+0),'FontSize',9);
text(xr(1)+0.08, S.loa_hi-0.12, sprintf('+1.96 SD = %+.2f', S.loa_hi),'FontSize',8.5,'Color',[0.3 0.3 0.3]);
text(xr(1)+0.08, S.loa_lo+0.12, sprintf('-1.96 SD = %+.2f', S.loa_lo),'FontSize',8.5,'Color',[0.3 0.3 0.3]);
set(gca,'FontSize',9); box on; hold off;
end

function S = same_subset_loocv(cfg, E_all, tps)
% Per-ROI LOOCV (loocv_err) restricted to the time points tps, and near (II+III) / remote (I+IV) pooled RMSE;
% *_all = the same pooling of the all-time-point held-out errors E_all (n_tp x 4 x method).
S = struct('name',{},'tps',{},'errors',{},'rmse_roi',{},'near',{},'remote',{},'near_all',{},'remote_all',{});
for q = 1:numel(cfg)
    Y = cfg{q}{1};  T = cfg{q}{2};
    E = loocv_err(Y(tps,:), T(tps,:));
    Ea = E_all(:,:,q);
    S(q).name = cfg{q}{3};  S(q).tps = tps;  S(q).errors = E;
    S(q).rmse_roi = sqrt(mean(E.^2,1));
    S(q).near = sqrt(mean(E(:,[2 3]).^2,'all'));    S(q).remote = sqrt(mean(E(:,[1 4]).^2,'all'));
    S(q).near_all = sqrt(mean(Ea(:,[2 3]).^2,'all'));  S(q).remote_all = sqrt(mean(Ea(:,[1 4]).^2,'all'));
end
end

function F = far_to_near(name, Y, Tp, t_q, ds, t_gre, alpha, scheme, n_cal)
% Far-to-near calibration of one 2DVSMSI variant with the GRE PRF signal as the far temperature reference.
%   Y, Tp  n_tp x 4 ROI-mean apparent-T2 deviation (ms) and probe temperature (degC; ROIs II, III only for scoring)
%   t_q    n_tp x 1 window midpoints of the 2DVSMSI scans; t_gre n_tp x 1 GRE window midpoints (scanner clock, s)
%   ds     n_tp x 4 GRE ROI-mean ds deviation; alpha literature PRF coefficient (ds per degC)
% For each test time point j the scheme fixes the calibration time points cal (= the "known" GRE scans):
%   fwd / fwd_nointerp: cal = 1..n_cal, j = n_cal+1..n_tp, anchor 1;   online_now: cal = 1..j, anchor 1;
%   online_past: cal = 1..j-1, anchor 1;   rev: cal = n_tp-n_cal+1..n_tp, j = 1..n_tp-n_cal, anchor n_tp.
% Far reference at ROIs I and IV: X = ds/alpha from the known GRE scans, linearly interpolated in time to t_q
% (linear extrapolation from the two nearest known points outside their range); fwd_nointerp uses ds/alpha at
% the same time point without interpolation.  m = common slope of Y on X over cal (ROIs I, IV, own intercepts).
% Near ROIs II and III: error = (Y(j)-Y(anchor))/m - (Tp(j)-Tp(anchor)).
n = size(Y,1);
switch scheme
    case {'fwd','fwd_nointerp'}, test = n_cal+1:n;  cal = arrayfun(@(j) 1:n_cal, test, 'UniformOutput', false);
    case 'online_now',           test = n_cal+1:n;  cal = arrayfun(@(j) 1:j,     test, 'UniformOutput', false);
    case 'online_past',          test = n_cal+1:n;  cal = arrayfun(@(j) 1:j-1,   test, 'UniformOutput', false);
    case 'rev',                  test = 1:n-n_cal;  cal = arrayfun(@(j) n-n_cal+1:n, test, 'UniformOutput', false);
    otherwise, error('stats:FarScheme','Unknown far-to-near scheme %s. Fix: use fwd, fwd_nointerp, online_now, online_past or rev.',scheme);
end
far = [1 4];  near = [2 3];
e = zeros(numel(test),2);  sl = zeros(1,numel(test));
for jj = 1:numel(test)
    c = cal{jj};  j = test(jj);
    X = zeros(n,4);
    for k = far
        if strcmp(scheme,'fwd_nointerp')
            X(:,k) = ds(:,k)/alpha;
        else
            X(:,k) = interp1(t_gre(c), ds(c,k)/alpha, t_q, 'linear', 'extrap');
        end
    end
    sl(jj) = anc(Y(c,far), X(c,far));
    if strcmp(scheme,'rev'), a = c(end); else, a = c(1); end
    e(jj,:) = (Y(j,near) - Y(a,near))/sl(jj) - (Tp(j,near) - Tp(a,near));
end
F.variant = name;  F.scheme = scheme;  F.test_tp = test;  F.anchor = a;  F.slopes = sl;  F.errors = e;
F.rmse = sqrt(mean(e.^2,1));  F.rmse_both = sqrt(mean(e.^2,'all'));
F.bias = mean(e,1);  F.maxe = max(abs(e),[],'all');
end

function P = roi_position(name, T2v, Mv, ctr, Tp, drs, dcs, ROD_R, ROD_C, PX)
% ROI II window (3x3) shifted by (dr, dc): ROI-mean apparent T2 over time, per-ROI slope (own intercept) on the
% probe temperature Tp, per-ROI LOOCV RMSE (loocv_err), distance to the rod surface (rod_dist, ROI II rule),
% mean 16 ms SOS magnitude relative to the unshifted window.  Voxels with T2 <= 0 or NaN are left out (count kept).
nd = numel(drs);  nc = numel(dcs);  n = size(T2v,3);
P.variant = name;  P.drow = drs;  P.dcol = dcs;
P.centre = zeros(1,2,nd,nc);  P.slope = zeros(nd,nc);  P.rmse = zeros(nd,nc);  P.dist = zeros(nd,nc);
P.mag = zeros(nd,nc);  P.n_invalid = zeros(nd,nc);  P.series = zeros(n,nd,nc);
for i = 1:nd
    for j = 1:nc
        c = ctr + [drs(i) dcs(j)];
        y = zeros(n,1);  mg = zeros(n,1);  bad = 0;
        for t = 1:n
            v = reshape(T2v(c(1)+(-1:1), c(2)+(-1:1), t),[],1);
            w = reshape(Mv(c(1)+(-1:1), c(2)+(-1:1), t),[],1);
            ok = isfinite(v) & v > 0;  bad = bad + nnz(~ok);
            y(t) = mean(v(ok));  mg(t) = mean(w(ok));
        end
        rf = roifit(y, Tp);
        P.centre(:,:,i,j) = c;  P.slope(i,j) = rf.slope;  P.rmse(i,j) = sqrt(mean(loocv_err(y, Tp).^2));
        P.dist(i,j) = rod_dist(c, ROD_R(1),ROD_R(2),ROD_C(1),ROD_C(2), PX, 2);
        P.mag(i,j) = mean(mg);  P.n_invalid(i,j) = bad;  P.series(:,i,j) = y;
    end
end
P.mag_rel = P.mag / P.mag(drs==0, dcs==0);
end

function G = pos_ranges(P, drs)
% Ranges over the windows shifted towards the rod end (drow < 0): slope and LOOCV per variant (rows = P),
% distance for drow = -1 and -2 and relative magnitude over both variants.
sh = drs < 0;
G.slope = zeros(numel(P),2);  G.rmse = zeros(numel(P),2);  d1 = [];  d2 = [];  mg = [];
for v = 1:numel(P)
    s = P(v).slope(sh,:);  r = P(v).rmse(sh,:);
    G.slope(v,:) = [min(s(:)) max(s(:))];  G.rmse(v,:) = [min(r(:)) max(r(:))];
    d1 = [d1, P(v).dist(drs==-1,:)];  d2 = [d2, P(v).dist(drs==-2,:)];  %#ok<AGROW>
    m = P(v).mag_rel(sh,:);  mg = [mg; m(:)];  %#ok<AGROW>
end
G.dist_one = [min(d1) max(d1)];  G.dist_two = [min(d2) max(d2)];  G.mag = [min(mg) max(mg)];
end

function P = twob_precision(B, m, NB, SEED)
% Single-measurement precision from the two back-to-back PFA blocks of the 2b scans, and block 1 as a separate
% PFA calibration series (common slope with cluster-bootstrap CI as S1b, per-ROI slopes, CV, per-ROI LOOCV as S2).
% B = load('rev_2b_blocks.mat') or [] (then every field is NaN and P.available = false).
P.available = ~isempty(B);  P.m = m;
if ~P.available
    P.r = nan(9,4);  P.dT_probe = nan(9,4);  P.dT_range = [NaN NaN];  P.n = NaN;
    P.sd_single = NaN;  P.sd_roi = nan(1,4);  P.mean_r = NaN;  P.sd_diff = NaN;
    P.b1_slopes = nan(1,4);  P.b1_r2 = nan(1,4);  P.b1_cv = NaN;  P.b1_common = NaN;  P.b1_ci = [NaN NaN];
    P.b1_errors = nan(9,4);  P.b1_rmse = nan(1,4);  P.b1_rmse_overall = NaN;
    return
end
r = (B.A1 - B.A2)/m - (B.T1 - B.T2);
P.r = r;  P.dT_probe = B.T1 - B.T2;  P.dT_range = [min(P.dT_probe(:)) max(P.dT_probe(:))];  P.n = numel(r);
P.sd_single = std(r(:))/sqrt(2);  P.sd_roi = std(r,0,1)/sqrt(2);  P.mean_r = mean(r(:));  P.sd_diff = std(r(:));
Y1 = B.A1 - mean(B.A1,1);                 % deviation from the time mean (the ROI intercepts absorb the reference)
[sl,mc] = common_fit(Y1, B.T1);  rf = roifit(Y1, B.T1);
[~,ci] = pooled_boot(Y1, B.T1, NB, SEED);
E = loocv_err(Y1, B.T1);
P.b1_slopes = sl;  P.b1_r2 = rf.r2;  P.b1_cv = 100*std(sl)/mean(sl);  P.b1_common = mc;  P.b1_ci = ci;
P.b1_errors = E;  P.b1_rmse = sqrt(mean(E.^2,1));  P.b1_rmse_overall = sqrt(mean(P.b1_rmse.^2));
end
