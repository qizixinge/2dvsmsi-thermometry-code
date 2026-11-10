% Purpose: generate manuscript-analysis inputs by reconstructing 18 raw 2DVSMSI files on the native 100x100 grid and importing nine GRE phase DICOM images.
%% reconstruct_all.m
% Single entry point for producing the raw inputs of the manuscript analysis.
% It replaces newexperiments_xulie.m (hand-written per-dataset blocks with
% manual ROI offsets) and rerun_native_100.m (2DVSMSI only).
%
% PART A  2DVSMSI: 9 time points x {PFA, Full} -> MIP and SOS T2 maps on the
%         NATIVE 100x100 grid.  No interpolation.  Writes
%         T2MAP_and_img4i_native100.mat.
% PART B  PRF/GRE : 9 phase DICOMs (acquisitions 1-9) -> phaseMAPr (128x128x9, radians), plus the
%         GRE magnitude image (series 0016 = acquisition 2; no other GRE magnitude exists).  Writes phaseMAPr_and_magmap.mat.
%
% ROI coordinates and probe temperatures live in ROI_and_tempr_resolved.mat
% and are NOT touched here; nothing in this file adds a manual offset.
%
% 2026-11-05: acquisition 1 added (time points 1-9 = acquisitions 1-9).
% Part A also stores fit_seconds (time point x [PFA Full] x [MIP SOS]), the
% RDNLS fitting time parsed from MAP's RDNLS_TIMING output.
%
% Runtime: ~10-20 min for Part A (POCS on 9 PFA datasets).  Part B is seconds.
% Set RUN_A / RUN_B below to skip a part whose output already exists.
%
% Created 2026-09-03 (merges rerun_native_100.m + the DICOM block of
% newexperiments_xulie.m, per plan section 5.3).

clear; close all;
project = project_paths();
cd(project.code);
addpath(project.mapvbvd);

RUN_A = true;    % 2DVSMSI reconstruction + T2 fitting
RUN_B = true;    % PRF phase DICOM loading

num = 9;                         % time points, cooling order = acquisitions 1-9
N   = 100;                       % native acquisition matrix
NECHO = 10;  NBINS = 12;         % echoes per bin, spectral bins

%% ========================================================================
%  PART A  2DVSMSI at native 100x100
%% ========================================================================
if RUN_A
    % PFA (pff = 0.6) and Full (pff = 1) file pairs, in cooling order
    % Acquisitions 1-9 of the session (acquisition 1 = MID00570/00572).
    files_pfa  = {'meas_MID00570_FID02828_pulseq.dat','meas_MID00576_FID02834_pulseq.dat','meas_MID00582_FID02840_pulseq.dat', ...
                  'meas_MID00588_FID02846_pulseq.dat','meas_MID00594_FID02852_pulseq.dat', ...
                  'meas_MID00600_FID02858_pulseq.dat','meas_MID00606_FID02864_pulseq.dat', ...
                  'meas_MID00612_FID02870_pulseq.dat','meas_MID00618_FID02876_pulseq.dat'};
    files_full = {'meas_MID00572_FID02830_pulseq.dat','meas_MID00578_FID02836_pulseq.dat','meas_MID00584_FID02842_pulseq.dat', ...
                  'meas_MID00590_FID02848_pulseq.dat','meas_MID00596_FID02854_pulseq.dat', ...
                  'meas_MID00602_FID02860_pulseq.dat','meas_MID00608_FID02866_pulseq.dat', ...
                  'meas_MID00614_FID02872_pulseq.dat','meas_MID00620_FID02878_pulseq.dat'};

    % 4th index: 1 = PFA, 3 = Full (index 2 unused, kept for compatibility
    % with the older 3-variant layout that downstream scripts expect)
    T2MAP    = zeros(N,N,num,3);
    T2MAPsos = zeros(N,N,num,3);
    img4i    = zeros(N,N,num,3);
    img4isos = zeros(N,N,num,3);
    fit_seconds = zeros(num,2,2);    % time point x [PFA Full] x [MIP SOS]

    fprintf('--- Part A: 2DVSMSI reconstruction at native %dx%d ---\n', N, N);
    for t = 1:num
        txt = evalc('[a,b,c,d] = MAP(NECHO,NBINS,0.6,fullfile(project.raw,files_pfa{t}),0)');
        fit_seconds(t,1,:) = parse_fit_seconds(txt);  fprintf('%s',txt);
        img4i(:,:,t,1)=a; img4isos(:,:,t,1)=b; T2MAP(:,:,t,1)=c; T2MAPsos(:,:,t,1)=d;
        txt = evalc('[a,b,c,d] = MAP(NECHO,NBINS,1,fullfile(project.raw,files_full{t}),0)');
        fit_seconds(t,2,:) = parse_fit_seconds(txt);  fprintf('%s',txt);
        img4i(:,:,t,3)=a; img4isos(:,:,t,3)=b; T2MAP(:,:,t,3)=c; T2MAPsos(:,:,t,3)=d;
        fprintf('  time point %d/%d done\n', t, num);
    end
    reconstruction_frame='centered-fft-v1';
    save(fullfile(project.data,'T2MAP_and_img4i_native100.mat'),'reconstruction_frame','T2MAP','T2MAPsos','img4i','img4isos','num','fit_seconds','files_pfa','files_full','-v7.3');
    fprintf('wrote T2MAP_and_img4i_native100.mat  (%dx%d x %d x 3)\n\n', N, N, num);
end

%% ========================================================================
%  PART B  PRF / GRE phase maps from DICOM
%% ========================================================================
if RUN_B
    % Phase DICOMs in cooling order (series 0015, 0017, 0019, ..., 0031; acquisitions 1-9).
    % Siemens stores phase as an unsigned integer; rescale to [-pi, pi).
    dicom_phase = { ...
      'GRE_series0015_phase.IMA', ...
      'GRE_series0017_phase.IMA', ...
      'GRE_series0019_phase.IMA', ...
      'GRE_series0021_phase.IMA', ...
      'GRE_series0023_phase.IMA', ...
      'GRE_series0025_phase.IMA', ...
      'GRE_series0027_phase.IMA', ...
      'GRE_series0029_phase.IMA', ...
      'GRE_series0031_phase.IMA'};
    % Magnitude of the GRE scan of acquisition 2 (phase series 0017).  No GRE
    % magnitude exists for the other acquisitions (series 0018, 0020, ... and the
    % acquisition-1 magnitude were not exported), so every GRE magnitude-based
    % quantity (display threshold, SNR, figure) uses this single image.
    dicom_mag = 'GRE_series0016_magnitude.IMA';

    fprintf('--- Part B: PRF phase maps ---\n');
    assert(numel(dicom_phase)==num,'reconstruct_all:Count','Expected %d phase DICOMs.',num);
    phaseMAPr = zeros(128,128,num);
    for t = 1:num
        ph = double(dicomread(fullfile(project.data,dicom_phase{t})));
        phaseMAPr(:,:,t) = ph/max(ph,[],'all')*2*pi - pi;
    end
    magmap = dicomread(fullfile(project.data,dicom_mag));
    save(fullfile(project.data,'phaseMAPr_and_magmap.mat'),'phaseMAPr','magmap');
    fprintf('wrote phaseMAPr_and_magmap.mat  (128x128 x %d)\n', num);
end

fprintf('\nNext: run statistics_analysis_final.m\n');


function s = parse_fit_seconds(txt)
% [MIP SOS] seconds from MAP's 'RDNLS_TIMING pff=.. MIP=.. s SOS=.. s' line.
tok = regexp(txt,'MIP=([\d.]+) s SOS=([\d.]+) s','tokens','once');
if isempty(tok), error('reconstruct_all:NoTiming','MAP printed no RDNLS_TIMING line. Fix: check MAP.m still prints it.'); end
s = str2double(tok);
end
