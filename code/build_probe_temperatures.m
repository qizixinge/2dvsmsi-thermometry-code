% Purpose: derive the probe temperature arrays (temp1, temp3, temprp1) from the fluoroptic log as the mean over each scan's acquisition window.
%% build_probe_temperatures.m
% Replaces the hand-recorded probe readings in ROI_and_tempr_resolved.mat
% with values computed from data/probe_log_20250425.csv (one sample per
% second, 0.1 degC resolution, channels CH1-CH4 = probes at ROI I-IV).
%
% Each entry is the mean of the log samples inside the scan's acquisition
% window.  Windows (scanner clock, 2025-04-25, acquisitions 1-9 = time points 1-9):
%   PFA / Full : first to last ADC time stamp of the raw .dat file
%                (PFA 29.7 s, Full 49.7 s)
%   GRE/PRF    : DICOM AcquisitionTime + TotalScanTimeSec (27 s)
% The thermometer clock ran 51 s behind the scanner clock (thermometer 90 s
% and scanner 39 s behind the same phone clock, both checked by the author),
% so every scanner time is shifted by -51 s before indexing the log.  The
% window covers data acquisition only (dummy scans before the first ADC are
% excluded for PFA/Full).  Seconds missing from the log are skipped.  The hand-recorded arrays are
% kept as temp1_manual / temp3_manual / temprp1_manual for provenance.
%
% 2026-11-05: acquisition 1 added back (9 time points).  Its GRE started ~5 min
% before its 2DVSMSI scans; that only mattered when one hand-recorded
% temperature was shared by all three scans.  Each scan now has its own probe
% temperature, so the acquisition is usable.  The per-time-point ROI series are
% identical across time points (ROI coordinates are fixed per ROI), so the 8-row
% series stored in ROI_and_tempr_resolved.mat are replicated to 9 rows.
% Also stores scan_start_s / scan_end_s (scanner clock, s after midnight; 9x3
% = GRE, PFA, Full), gre_to_pfa_gap_s, and the within-window probe statistics
% quoted in the paper (win_range_max, win_dev_max).
%
% Run once after changing the time-point selection; statistics_analysis_final.m
% then reads the updated arrays.

clear;
project = project_paths();

fid = fopen(fullfile(project.data,'probe_log_20250425.csv'),'r');
fgetl(fid);                                          % header (GBK text)
C = textscan(fid,'%s %s %f %f %f %f %*[^\n]','Delimiter',',','EmptyValue',NaN);
fclose(fid);
tsec = cellfun(@(s) [3600 60 1]*sscanf(s,'%d:%d:%d'), C{2});
vals = [C{3} C{4} C{5} C{6}];
ok   = all(isfinite(vals),2) & any(vals>0,2);
tsec = tsec(ok);  vals = vals(ok,:);

hms = @(s) [3600 60 1]*sscanf(s,'%d:%d:%f');
% acquisitions 1-9, cooling order
gre_start = {'17:03:10.285','17:13:43.26','17:20:26.28','17:26:39.25','17:34:52.31', ...
             '17:42:56.28','17:49:18.32','17:56:44.30','18:05:09.25'};
gre_dur   = 27;
pfa_win   = {'17:08:12.32','17:08:42.00'; '17:14:16.40','17:14:46.08'; '17:21:00.65','17:21:30.33'; ...
             '17:27:12.46','17:27:42.14'; '17:35:27.39','17:35:57.07'; ...
             '17:43:29.43','17:43:59.11'; '17:49:53.43','17:50:23.11'; ...
             '17:57:20.43','17:57:50.11'; '18:05:42.38','18:06:12.06'};
full_win  = {'17:09:55.46','17:10:45.14'; '17:16:03.43','17:16:53.11'; '17:22:51.43','17:23:41.11'; ...
             '17:29:02.43','17:29:52.11'; '17:37:13.43','17:38:03.11'; ...
             '17:45:16.43','17:46:06.11'; '17:51:45.42','17:52:35.10'; ...
             '17:59:07.39','17:59:57.07'; '18:07:26.43','18:08:16.11'};

clock_offset = -51;                                  % thermometer = scanner - 51 s
winmean = @(t0,t1) mean(vals(tsec>=ceil(t0+clock_offset) & tsec<=floor(t1+clock_offset),:),1);
num = 9;
assert(numel(gre_start)==num && size(pfa_win,1)==num && size(full_win,1)==num, ...
    'build_probe_temperatures:Count','Expected %d scan windows per sequence. Fix: complete the gre_start/pfa_win/full_win lists.',num);
temprp1 = zeros(num,4);  temp1 = zeros(num,4);  temp3 = zeros(num,4);
scan_start_s = zeros(num,3);  scan_end_s = zeros(num,3);
win_range_max = 0;  win_dev_max = 0;
for t = 1:num
    g = hms(gre_start{t});
    w = [g g+gre_dur; hms(pfa_win{t,1}) hms(pfa_win{t,2}); hms(full_win{t,1}) hms(full_win{t,2})];
    scan_start_s(t,:) = w(:,1)';  scan_end_s(t,:) = w(:,2)';
    mm = zeros(3,4);
    for q = 1:3
        in = tsec>=ceil(w(q,1)+clock_offset) & tsec<=floor(w(q,2)+clock_offset);
        v = vals(in,:);  mm(q,:) = mean(v,1);
        win_range_max = max(win_range_max, max(max(v,[],1)-min(v,[],1)));
        win_dev_max   = max(win_dev_max,   max(abs(v-mm(q,:)),[],'all'));
    end
    temprp1(t,:) = mm(1,:);  temp1(t,:) = mm(2,:);  temp3(t,:) = mm(3,:);
end
gre_to_pfa_gap_s = scan_start_s(:,2) - scan_start_s(:,1);   % GRE start -> first PFA ADC
fprintf('GRE start -> first PFA ADC gap (s): %s\n', mat2str(gre_to_pfa_gap_s',4));
fprintf('within-window probe range max %.2f degC, max |sample - window mean| %.2f degC\n', win_range_max, win_dev_max);

f = fullfile(project.data,'ROI_and_tempr_resolved.mat');
R = load(f);
if ~isfield(R,'temp1_manual')          % keep the original hand-recorded values once
    R.temp1_manual = R.temp1;  R.temp3_manual = R.temp3;  R.temprp1_manual = R.temprp1;
end
R.temp1 = temp1;  R.temp3 = temp3;  R.temprp1 = temprp1;
R.num = num;
R.scan_start_s = scan_start_s;  R.scan_end_s = scan_end_s;  R.gre_to_pfa_gap_s = gre_to_pfa_gap_s;
R.win_range_max = win_range_max;  R.win_dev_max = win_dev_max;
% ROI coordinates are constant over time points: replicate the series to num rows
roi_fix = {{'ROI_2DVSMSI_100_series','ROI_2DVSMSI_100'},{'ROI_PRF_128_series','ROI_PRF_128'}};
for q = 1:numel(roi_fix)
    ser = R.(roi_fix{q}{1});  ref = double(R.(roi_fix{q}{2}));
    for k = 1:numel(ser)
        rows = double(ser{k});
        assert(all(rows == rows(1,:),'all') && all(rows(1,:) == ref(k,:)), 'build_probe_temperatures:RoiSeries', ...
            'ROI series %s{%d} is not constant over time points or differs from %s. Fix: supply the per-time-point coordinates.', roi_fix{q}{1}, k, roi_fix{q}{2});
        ser{k} = repmat(ref(k,:), num, 1);
    end
    R.(roi_fix{q}{1}) = ser;
end
R.temperature_source = ['probe_log_20250425.csv, mean over each scan acquisition window; ' ...
    'thermometer clock = scanner clock - 51 s; acquisitions 1-9 of 2025-04-25 = time points 1-9 (acquisition 1 re-included 2026-11-05: each scan has its own probe temperature). ' ...
    'temp*_manual = earlier hand-recorded readings of acquisitions 1-8 (8 rows), superseded.'];
save(f,'-struct','R');
fprintf('updated %s\n', f);
disp('temprp1 (GRE/PRF)'); disp(temprp1);
disp('temp1 (PFA)');       disp(temp1);
disp('temp3 (Full)');      disp(temp3);
