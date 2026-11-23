% Readout and phase-encoding FOV of the scanned 2DVSMSI sequences, from the
% k-space trajectory that Pulseq computes for the .seq files (calculateKspacePP).
% Needs the Pulseq MATLAB toolbox and the two scanned .seq files. They are looked up in
% this order: VSMSI_SEQ_DIR (if set), the repository (sequence/ for the .seq files,
% third_party/pulseq/matlab for Pulseq), then the parent folder of the repository.
here = fileparts(mfilename('fullpath')); root = fileparts(here);
addpath(root); project = project_paths();
files = {'writemy4TE1b_100_SHARE_full.seq','writemy4TE1b_100_SHARE_pocs06.seq'};
cand = {getenv('VSMSI_SEQ_DIR'), fullfile(project.root,'sequence'), fileparts(project.root)};
cand = cand(~cellfun(@isempty,cand));
seqdir = '';
for c_ = 1:numel(cand)
    if all(cellfun(@(f) isfile(fullfile(cand{c_},f)), files)), seqdir = cand{c_}; break; end
end
if isempty(seqdir)
    error('rev:seq',['%s and %s not found in: %s. Fix: copy them into %s or set VSMSI_SEQ_DIR to the folder ' ...
        'that holds them.'], files{1}, files{2}, strjoin(cand,' | '), fullfile(project.root,'sequence'));
end
pmc = {fullfile(project.root,'third_party','pulseq','matlab'), fullfile(seqdir,'pulseq-master','pulseq-master','matlab')};
pm = '';
for c_ = 1:numel(pmc)
    if isfile(fullfile(pmc{c_},'+mr','@Sequence','Sequence.m')), pm = pmc{c_}; break; end
end
if isempty(pm)
    error('rev:pulseq',['Pulseq toolbox not found in: %s. Fix: copy the matlab folder of Pulseq (v1.4) to ' ...
        '%s or set VSMSI_SEQ_DIR to a folder that contains pulseq-master.'], strjoin(pmc,' | '), pmc{1});
end
addpath(pm);
fprintf('Pulseq: %s\n.seq files: %s\n', pm, seqdir);
% writemy4TE1b_100_SHARE_pocs.seq (5 shots, ky -5..49) is an earlier design;
% pocs06 matches the scanned PFA protocol (6 shots, ky -10..49, 35 s).
FOV = struct();
for i = 1:2
    f = fullfile(seqdir,files{i});
    if ~exist(f,'file'), error('rev:seq','%s not found. Set VSMSI_SEQ_DIR.',f); end
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
if ~isfolder(outdir), mkdir(outdir); end
save(fullfile(outdir,'readout_fov.mat'),'FOV');
