% Pseudo-multiple-replica noise simulation with the actual sampling and reconstruction.
% Complex Gaussian noise with the measured channel covariance is added to the
% acquired k-space of time point PR_TP (default 3); each replica is reconstructed with
% rev_recon_channels (k-space ordering of dataSOS/dataprocessing, POCS for PFA),
% combined by SOS and MIP, and fitted with the RDNLS grid of RDNLS.m.
% Inputs (optional, set before running): PR_NR (replicas per chunk), PR_CHUNK,
% PR_VARIANTS (subset of {'PFA','Full'}), PR_TP (time point index).
% Output: VSMSI_REV_DIR/pseudoreplica/<variant>_tp<k>_chunk<j>.mat
here = fileparts(mfilename('fullpath')); root = fileparts(here);
addpath(root); addpath(here);
project = project_paths(); addpath(project.mapvbvd);
outdir = getenv('VSMSI_REV_DIR'); if isempty(outdir), outdir = project.supporting; end
prdir = fullfile(outdir,'pseudoreplica'); if ~isfolder(prdir), mkdir(prdir); end
if ~exist('PR_NR','var'), PR_NR = 100; end
if ~exist('PR_CHUNK','var'), PR_CHUNK = 1; end
if ~exist('PR_VARIANTS','var'), PR_VARIANTS = {'PFA','Full'}; end
if ~exist('PR_TP','var'), PR_TP = 3; end     % time point of the illustrated images (TP_CMP)

% Acquisitions 1-9 of the session (= time points 1-9 of the paper).
files_pfa  = {'meas_MID00570_FID02828_pulseq.dat','meas_MID00576_FID02834_pulseq.dat','meas_MID00582_FID02840_pulseq.dat', ...
              'meas_MID00588_FID02846_pulseq.dat','meas_MID00594_FID02852_pulseq.dat', ...
              'meas_MID00600_FID02858_pulseq.dat','meas_MID00606_FID02864_pulseq.dat', ...
              'meas_MID00612_FID02870_pulseq.dat','meas_MID00618_FID02876_pulseq.dat'};
files_full = {'meas_MID00572_FID02830_pulseq.dat','meas_MID00578_FID02836_pulseq.dat','meas_MID00584_FID02842_pulseq.dat', ...
              'meas_MID00590_FID02848_pulseq.dat','meas_MID00596_FID02854_pulseq.dat', ...
              'meas_MID00602_FID02860_pulseq.dat','meas_MID00608_FID02866_pulseq.dat', ...
              'meas_MID00614_FID02872_pulseq.dat','meas_MID00620_FID02878_pulseq.dat'};

% ---- noise model: channel covariance from noise-only spectral bins of the
%      fully sampled acquisition, central readout rows (flat receiver passband)
nm_file = fullfile(prdir,'noise_model_bins.mat');
if ~isfile(nm_file)
    tw = mapVBVD(fullfile(project.raw,files_full{PR_TP})); raw = tw{end}.image(); clear tw
    ims = rev_recon_channels(raw,10,1,12,1);
    c = ims(:,:,:,:,1);                                  % TE index 1; the 4 TE images share peripheral lines
    E = squeeze(sum(abs(c).^2,[1 2 3]));
    [~,o] = sort(E); noise_bins = sort(o(1:3))';
    rows_flat = 25:75;
    x = reshape(c(:,rows_flat,:,noise_bins),size(c,1),[]);
    Psi_img = (x*x')/size(x,2);
    Psi_k = Psi_img*numel(c(1,:,:,1));                 % ifft2 scaling: var_img = var_k/(Nx*Ny)
    prof = squeeze(mean(abs(c(:,:,:,noise_bins)).^2,[1 3 4]));
    save(nm_file,'Psi_img','Psi_k','noise_bins','rows_flat','prof','E');
    fprintf('noise model: bins %s (relative energy %s), rows %d-%d, %d samples\n', mat2str(noise_bins), ...
        mat2str(round(E(noise_bins)'/min(E),2)), rows_flat(1), rows_flat(end), size(x,2));
end
load(nm_file,'Psi_k');
L = chol((Psi_k+Psi_k')/2,'lower');

for v = 1:numel(PR_VARIANTS)
    vn = PR_VARIANTS{v};
    if strcmp(vn,'PFA'), fn = files_pfa{PR_TP}; pff = 0.6; else, fn = files_full{PR_TP}; pff = 1; end
    tw = mapVBVD(fullfile(project.raw,fn)); raw0 = tw{end}.image(); clear tw
    [nx,nc,nl] = size(raw0);
    seed = 20261012 + 1000*strcmp(vn,'Full') + 100*PR_TP + PR_CHUNK;
    rng(seed,'twister');
    sos3 = zeros(100,100,3,PR_NR,'single'); mip3 = sos3;
    T2s = zeros(100,100,PR_NR,'single'); T2m = T2s;
    t0 = tic;
    for r = 1:PR_NR
        w = L*(randn(nc,nx*nl)+1i*randn(nc,nx*nl))/sqrt(2);
        w = permute(reshape(w,nc,nx,nl),[2 1 3]);
        ims = rev_recon_channels(raw0+w,10,1,12,pff);
        mag = squeeze(sqrt(sum(abs(ims).^2,1)));
        s = squeeze(sqrt(sum(mag(:,:,:,2:4).^2,3)));
        m = squeeze(max(mag(:,:,:,2:4),[],3));
        sos3(:,:,:,r) = s; mip3(:,:,:,r) = m;
        T2s(:,:,r) = reshape(rev_t2fit(reshape(s,[],3)),100,100);
        T2m(:,:,r) = reshape(rev_t2fit(reshape(m,[],3)),100,100);
    end
    fprintf('%s tp%d chunk %d: %d replicas in %.0f s (seed %d)\n', vn, PR_TP, PR_CHUNK, PR_NR, toc(t0), seed);
    save(fullfile(prdir,sprintf('%s_tp%d_chunk%d.mat',vn,PR_TP,PR_CHUNK)),'sos3','mip3','T2s','T2m','seed','pff','fn','-v7.3');
end
