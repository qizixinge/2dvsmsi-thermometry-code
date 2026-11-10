% Readout and phase-encoding FOV of the scanned 2DVSMSI sequences, from the
% k-space trajectory that Pulseq computes for the .seq files (calculateKspacePP).
% Needs the Pulseq MATLAB toolbox (set VSMSI_PULSEQ_DIR to its matlab/ folder, i.e. the folder
% that contains +mr) and the two .seq files written by the scripts in sequence/ (run
% sequence/writemy4TE1b_100_SHARE_full.m and sequence/write4TE1b_100_SHARE_pocs.m first;
% set VSMSI_SEQ_DIR to use another folder).
here = fileparts(mfilename('fullpath')); root = fileparts(here);
addpath(root); project = project_paths();
seqdir = getenv('VSMSI_SEQ_DIR'); if isempty(seqdir), seqdir = fullfile(project.root,'sequence'); end
pm = getenv('VSMSI_PULSEQ_DIR');
if isempty(pm) || ~exist(fullfile(pm,'+mr','@Sequence','Sequence.m'),'file')
    error('rev:pulseq','Pulseq toolbox not found. Fix: set VSMSI_PULSEQ_DIR to the matlab/ folder of the Pulseq toolbox (it contains +mr).');
end
addpath(pm);
% writemy4TE1b_100_SHARE_pocs.seq (5 shots, ky -5..49) is an earlier design;
% pocs06 matches the scanned PFA protocol (6 shots, ky -10..49, 35 s).
files = {'writemy4TE1b_100_SHARE_full.seq','writemy4TE1b_100_SHARE_pocs06.seq'};
FOV = struct();
for i = 1:2
    f = fullfile(seqdir,files{i});
    if ~exist(f,'file'), error('rev:seq','%s not found. Fix: run the sequence/ writer scripts or set VSMSI_SEQ_DIR.',f); end
    seq = mr.Sequence(); seq.read(f);
    ktraj_adc = seq.calculateKspacePP();
    kx = ktraj_adc(1,:); ky = ktraj_adc(2,:);
    nadc = 100; K = reshape(kx,nadc,[]); Y = reshape(ky,nadc,[]);
    dkx = diff(K); dkx = dkx(:);
    kyl = unique(round(mean(Y),6)); dky = min(diff(kyl));
    FOV(i).file = files{i};
    FOV(i).n_lines = size(K,2);
    FOV(i).dkx = median(abs(dkx)); FOV(i).dkx_range = [min(abs(dkx)) max(abs(dkx))];
    FOV(i).kx_span = median(max(K)-min(K));
    FOV(i).dky = dky; FOV(i).ky_index = round(kyl/dky);
    FOV(i).fov_ro_mm = 1e3/FOV(i).dkx; FOV(i).fov_pe_mm = 1e3/dky;
    fprintf('%s: %d ADC lines; readout dk %.5f 1/m (range %.5f-%.5f) -> FOV_RO %.2f mm, pixel %.4f mm; PE dk %.5f 1/m -> FOV_PE %.2f mm, pixel %.4f mm; ky index %d..%d (%d unique)\n', ...
        files{i}, FOV(i).n_lines, FOV(i).dkx, FOV(i).dkx_range, FOV(i).fov_ro_mm, FOV(i).fov_ro_mm/100, ...
        dky, FOV(i).fov_pe_mm, FOV(i).fov_pe_mm/100, min(FOV(i).ky_index), max(FOV(i).ky_index), numel(kyl));
end
fprintf('analytic: FOV_RO = 256 mm x readoutTime/samplingTime = 256 x 0.54/0.52 = %.2f mm\n', 256*0.54/0.52);
outdir = getenv('VSMSI_REV_DIR'); if isempty(outdir), outdir = project.supporting; end
save(fullfile(outdir,'readout_fov.mat'),'FOV');
