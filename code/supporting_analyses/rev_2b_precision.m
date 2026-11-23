% Purpose: reconstruct the two back-to-back PFA blocks of the nine 2025-04-25 "2b" scans and save their ROI means and probe temperatures for the single-measurement precision.
%% rev_2b_precision.m
% The 2025-04-25 "2b" scans (MID00571-MID00619, every 6th MID; sequence writemy4TE2b_100_SHARE_pocs06.seq)
% acquire two complete PFA data sets back to back in one file:
%   block 1 = imaging ADCs   1..936  (one dummy shot without ADC, then the same shots as the 1b PFA scan)
%   block 2 = imaging ADCs 937..1872 (the same 60 phase-encoding lines again, no dummy shot)
% No k-space line is shared between the blocks, so the two blocks are two independent PFA measurements a
% few seconds apart.  statistics_analysis_final.m uses their difference for the single-measurement
% precision and block 1 as a second PFA calibration series.
%
% Each block is reconstructed with the paper's pipeline (MAP.m): mapVBVD reader (image.unsorted() gives the
% imaging ADCs in acquisition order, which equals image() for these unlabeled Pulseq scans; checked below on
% the 1b file when it is available), dataSOS with partial-Fourier factor 0.6 (POCS), first echo (8 ms)
% removed, RDNLS on the 16/24/32 ms echoes (tesp 8 ms, T2approx 60).  ROI means: 3x3 voxels at the paper's
% PFA ROI centres (analysis_geometry pfa_roi).  Probe temperature of each block: mean of the 1 Hz log
% (data/probe_log_20250425.csv) over the block's first..last ADC time stamp, thermometer clock = scanner
% clock - 51 s, the same rule as build_probe_temperatures.m.
%
% Input : raw files meas_MID00571_*.dat ... meas_MID00619_*.dat in VSMSI_RAW_DIR (default data/).
%         They are not part of the repository (restricted raw data).
% Output: data/rev_2b_blocks.mat (small derived file, committed; statistics_analysis_final.m reads it).
% WARNING: this script OVERWRITES data/rev_2b_blocks.mat.
% Run time: about 1 min per block (18 blocks).
here = fileparts(mfilename('fullpath')); root = fileparts(here);
addpath(root);
project = project_paths(); addpath(project.mapvbvd);

MIDS    = 571:6:619;            % acquisitions 1-9 of 2025-04-25 (time points 1-9)
NPER    = 936;                  % imaging ADCs per block (13 echoes x 12 bins x 6 shots)
NECHO   = 10;  NBINS = 12;  PFF = 0.6;
CLOCK_OFFSET = -51;             % thermometer = scanner - 51 s (build_probe_temperatures.m)
SEQNAME = 'writemy4TE2b_100_SHARE_pocs06.seq';
n_tp = numel(MIDS);

% ---- raw files ------------------------------------------------------------------------------------
roots = unique({project.raw, project.data}, 'stable');
files = cell(1,n_tp);
for t = 1:n_tp
    files{t} = find_raw(MIDS(t), roots);
    if isempty(files{t})
        error('rev2b:MissingRaw', ['Raw file meas_MID%05d_*.dat (2025-04-25 "2b" PFA scan, %s) was not found in %s. ' ...
            'Fix: put the nine raw files meas_MID00571_*.dat ... meas_MID00619_*.dat (every 6th MID) into one folder and ' ...
            'set the environment variable VSMSI_RAW_DIR to it before running this script, e.g. ' ...
            'setenv(''VSMSI_RAW_DIR'',''D:\\path\\to\\raw''), or copy them into data/. ' ...
            'statistics_analysis_final.m does not need them: it reads the committed data/rev_2b_blocks.mat.'], ...
            MIDS(t), SEQNAME, strjoin(roots,' and '));
    end
end

% ---- geometry and probe log ------------------------------------------------------------------------
S1f = load(fullfile(project.data,'T2MAP_and_img4i_native100.mat'),'T2MAPsos','reconstruction_frame');
R   = load(fullfile(project.data,'ROI_and_tempr_resolved.mat'));
g   = analysis_geometry(R, S1f.reconstruction_frame);
roi = g.pfa_roi;                                   % 4x2 [row col], native grid, PFA frame
fid = fopen(fullfile(project.data,'probe_log_20250425.csv'),'r');
fgetl(fid);                                        % header (GBK text)
C = textscan(fid,'%s %s %f %f %f %f %*[^\n]','Delimiter',',','EmptyValue',NaN);
fclose(fid);
tsec = cellfun(@(s) [3600 60 1]*sscanf(s,'%d:%d:%d'), C{2});
vals = [C{3} C{4} C{5} C{6}];
ok   = all(isfinite(vals),2) & any(vals>0,2);
tsec = tsec(ok);  vals = vals(ok,:);
winmean = @(t0,t1) mean(vals(tsec>=ceil(t0+CLOCK_OFFSET) & tsec<=floor(t1+CLOCK_OFFSET),:),1);
winn    = @(t0,t1) nnz(tsec>=ceil(t0+CLOCK_OFFSET) & tsec<=floor(t1+CLOCK_OFFSET));

% ---- reader/pipeline check on the 1b PFA file of time point 1 (if present) -------------------------
% unsorted() must equal image() and the block pipeline must reproduce the committed PFA-SOS map and the
% committed PFA window start of time point 1.
f1b = find_raw(570, roots);
check_1b = struct('done',false,'max_raw_diff',NaN,'max_t2_diff',NaN,'start_diff_s',NaN);
if ~isempty(f1b)
    tw = mapVBVD(f1b);  if iscell(tw), tw = tw{end}; end
    raw_u = tw.image.unsorted();  raw_s = tw.image();
    check_1b.max_raw_diff = max(abs(raw_u(:) - raw_s(:)));
    [T2s,~] = recon_block(raw_u, NECHO, NBINS, PFF);
    check_1b.max_t2_diff = max(abs(T2s - double(S1f.T2MAPsos(:,:,1,1))),[],'all');
    ts1 = double(tw.image.timestamp(:))*0.0025;
    check_1b.start_diff_s = ts1(1) - R.scan_start_s(1,2);
    check_1b.done = true;
    fprintf('1b check (%s): max|unsorted - image| %.3g, max|T2 - committed PFA-SOS TP1| %.3g ms, window start - scan_start_s %.3g s\n', ...
        f1b, check_1b.max_raw_diff, check_1b.max_t2_diff, check_1b.start_diff_s);
    assert(check_1b.max_raw_diff == 0 && check_1b.max_t2_diff < 1e-9 && abs(check_1b.start_diff_s) < 1e-6, 'rev2b:ReaderCheck', ...
        'The unsorted reader or the block pipeline does not reproduce the committed PFA time point 1. Fix: check mapVBVD and recon_block against MAP.m before using the 2b results.');
    clear tw raw_u raw_s
else
    fprintf('1b check skipped: meas_MID00570_*.dat not found in %s\n', strjoin(roots,' and '));
end

% ---- reconstruct both blocks of every 2b scan ---------------------------------------------------------
A1 = zeros(n_tp,4);  A2 = zeros(n_tp,4);           % ROI-mean apparent T2 (ms), block 1 / block 2
T1 = zeros(n_tp,4);  T2 = zeros(n_tp,4);           % probe temperature over each block window (degC)
t_start = zeros(n_tp,2);  t_end = zeros(n_tp,2);   % first / last ADC time stamp, scanner clock (s after midnight)
n_log = zeros(n_tp,2);  seqnames = cell(1,n_tp);  nadc = zeros(1,n_tp);
for t = 1:n_tp
    tic_ = tic;
    tw = mapVBVD(files{t});  if iscell(tw), tw = tw{end}; end
    seqnames{t} = header_seqname(tw);
    if ~isempty(seqnames{t}) && ~strcmp(seqnames{t}, SEQNAME)
        error('rev2b:Sequence','%s was acquired with %s, not %s. Fix: check the MID list.', files{t}, seqnames{t}, SEQNAME);
    end
    raw = tw.image.unsorted();                     % [Col Cha NAcq], acquisition order
    ts  = double(tw.image.timestamp(:))*0.0025;    % 2.5 ms ticks -> s
    nadc(t) = size(raw,3);
    assert(nadc(t) == 2*NPER && numel(ts) == nadc(t), 'rev2b:Nadc', ...
        '%s has %d imaging ADCs (expected %d = two PFA blocks). Fix: check that the file is a 2b pocs06 scan.', files{t}, nadc(t), 2*NPER);
    for b = 1:2
        idx = (b-1)*NPER + (1:NPER);
        T2s = recon_block(raw(:,:,idx), NECHO, NBINS, PFF);
        A = zeros(1,4);
        for k = 1:4
            A(k) = mean(T2s(roi(k,1)+(-1:1), roi(k,2)+(-1:1)),'all');
        end
        t_start(t,b) = ts(idx(1));  t_end(t,b) = ts(idx(end));
        Tb = winmean(t_start(t,b), t_end(t,b));  n_log(t,b) = winn(t_start(t,b), t_end(t,b));
        assert(n_log(t,b) >= 20, 'rev2b:ProbeWindow', 'Only %d probe-log samples in the window of MID %d block %d. Fix: check the time stamps and the clock offset.', n_log(t,b), MIDS(t), b);
        if b == 1, A1(t,:) = A; T1(t,:) = Tb; else, A2(t,:) = A; T2(t,:) = Tb; end
    end
    fprintf('MID %d: block 1 %s-%s, block 2 %s-%s; T2app b1 %s b2 %s ms (%.0f s)\n', MIDS(t), hms(t_start(t,1)), hms(t_end(t,1)), ...
        hms(t_start(t,2)), hms(t_end(t,2)), sprintf('%.3f ',A1(t,:)), sprintf('%.3f ',A2(t,:)), toc(tic_));
    clear tw raw
end
assert(all(t_start(:,2) > t_end(:,1)), 'rev2b:BlockOrder', 'Block 2 does not start after block 1. Fix: check the ADC split.');

[~, fn, fe] = cellfun(@fileparts, files, 'UniformOutput', false);
raw_files = strcat(fn, fe);
mids = MIDS;  roi_pfa = roi;  n_per_block = NPER;  clock_offset_s = CLOCK_OFFSET;
reconstruction = sprintf(['mapVBVD image.unsorted(), imaging ADCs split 1..%d / %d..%d; dataSOS(raw,%d,1,%d,%.1f) (POCS); ' ...
    'first echo removed; RDNLS(.,[1 2 3],8,60) as MAP.m; ROI mean = 3x3 voxels at analysis_geometry pfa_roi'], NPER, NPER+1, 2*NPER, NECHO, NBINS, PFF);
probe_rule = sprintf(['mean of probe_log_20250425.csv samples (rows with four finite values, at least one > 0) with ' ...
    'ceil(t_start%+d) <= t_log <= floor(t_end%+d); t_start/t_end = first/last ADC time stamp of the block (scanner clock)'], CLOCK_OFFSET, CLOCK_OFFSET);
provenance = sprintf('code/supporting_analyses/rev_2b_precision.m, %s; sequence %s; MATLAB %s', ...
    datestr(now,'yyyy-mm-dd'), SEQNAME, version);
save(fullfile(project.data,'rev_2b_blocks.mat'), 'A1','A2','T1','T2','t_start','t_end','n_log','mids','raw_files','seqnames', ...
    'nadc','n_per_block','roi_pfa','clock_offset_s','reconstruction','probe_rule','provenance','check_1b');
fprintf('wrote %s\n', fullfile(project.data,'rev_2b_blocks.mat'));

function f = find_raw(mid, roots)
f = '';
for r = 1:numel(roots)
    d = dir(fullfile(roots{r}, sprintf('meas_MID%05d_*.dat', mid)));
    if ~isempty(d), f = fullfile(d(1).folder, d(1).name); return; end
end
end

function [T2sos, img] = recon_block(raw, necho, bins, pff)
% SOS branch of MAP.m on one block of imaging ADCs: dataSOS, drop the 8 ms echo, RDNLS on 16/24/32 ms.
sos = dataSOS(raw, necho, 1, bins, pff);
sos(:,:,1) = [];
T2sos = RDNLS(sos, [1 2 3], 8, 60);
img = sos(:,:,1);
end

function s = header_seqname(tw)
% Pulseq sequence file name stored by the scanner interpreter (empty if the header field is absent).
s = '';
try
    s = char(strtrim(strrep(string(tw.hdr.MeasYaps.sWipMemBlock.tFree),'"','')));   % char, so that the .mat stays readable outside MATLAB
catch
end
end

function s = hms(t)
s = sprintf('%02d:%02d:%05.2f', floor(t/3600), floor(mod(t,3600)/60), mod(t,60));
end
