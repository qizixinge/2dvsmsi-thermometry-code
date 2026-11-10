% Check that rev_recon_channels + rev_t2fit reproduce the saved native-grid maps.
here = fileparts(mfilename('fullpath')); root = fileparts(here);
addpath(root); addpath(here);
project = project_paths(); addpath(project.mapvbvd);
S1f = load(fullfile(project.data,'T2MAP_and_img4i_native100.mat'),'T2MAP','T2MAPsos','img4i','img4isos');
% time points 1 (acquisition 1, newly added) and 3 (illustrated time point)
files = {'meas_MID00570_FID02828_pulseq.dat', 0.6, 1, 1; 'meas_MID00572_FID02830_pulseq.dat', 1, 3, 1; ...
         'meas_MID00582_FID02840_pulseq.dat', 0.6, 1, 3; 'meas_MID00584_FID02842_pulseq.dat', 1, 3, 3};
for f = 1:size(files,1)
    tw = mapVBVD(fullfile(project.raw,files{f,1})); raw = tw{end}.image(); clear tw
    tic; ims = rev_recon_channels(raw,10,1,12,files{f,2}); trec = toc;
    mag = squeeze(sqrt(sum(abs(ims).^2,1)));          % Nx Nyo bins 4
    sos = squeeze(sqrt(sum(mag.^2,3)));               % Nx Nyo 4
    mip = squeeze(max(mag,[],3));
    v = files{f,3};  tpi = files{f,4};
    d_sos = max(abs(sos(:,:,2) - S1f.img4isos(:,:,tpi,v)),[],'all')/max(sos(:,:,2),[],'all');
    d_mip = max(abs(mip(:,:,2) - S1f.img4i(:,:,tpi,v)),[],'all')/max(mip(:,:,2),[],'all');
    tic; t2s = reshape(rev_t2fit(reshape(sos(:,:,2:4),[],3)),100,100); tfit = toc;
    t2m = reshape(rev_t2fit(reshape(mip(:,:,2:4),[],3)),100,100);
    ref_s = S1f.T2MAPsos(:,:,tpi,v); ref_m = S1f.T2MAP(:,:,tpi,v);
    ok_s = ref_s>0; ok_m = ref_m>0;
    fprintf('tp %d pff %.1f: recon %.1f s, fit %.1f s | rel max diff TE16 SOS %.2e MIP %.2e | T2 SOS max|d| %.3g ms (n=%d) MIP %.3g ms (n=%d)\n', ...
        tpi, files{f,2}, trec, tfit, d_sos, d_mip, max(abs(t2s(ok_s)-ref_s(ok_s))), nnz(ok_s), max(abs(t2m(ok_m)-ref_m(ok_m))), nnz(ok_m));
end
