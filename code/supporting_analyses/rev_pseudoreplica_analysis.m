% Summarise the pseudo-replica noise: spatial correlation of image and apparent-T2
% noise, cross-TE correlation, the noise SD of the 3x3 ROI mean, and
% replica-based SNR. Also evaluates the
% linear view-sharing model for the fully sampled acquisition.
here = fileparts(mfilename('fullpath')); root = fileparts(here);
addpath(root); addpath(here);
if ~exist('PR_TP','var'), PR_TP = 3; end     % time point of the replicas (rev_pseudoreplica.m)
if ~exist('PR_NR_MAX','var'), PR_NR_MAX = 200; end   % replicas used per variant (chunk order, first PR_NR_MAX)
project = project_paths();
outdir = getenv('VSMSI_REV_DIR'); if isempty(outdir), outdir = project.supporting; end
prdir = fullfile(outdir,'pseudoreplica');
S1f = load(fullfile(project.data,'T2MAP_and_img4i_native100.mat'),'T2MAP','T2MAPsos','img4isos','img4i','reconstruction_frame');
R = load(fullfile(project.data,'ROI_and_tempr_resolved.mat'));
S3f = load(fullfile(project.data,'YX.mat'),'yangaib');
ST = load(fullfile(project.results,'statistics_results_final.mat'),'POOL');
g = analysis_geometry(R,S1f.reconstruction_frame);
m_pool = [ST.POOL(1).pooled ST.POOL(2).pooled];          % PFA-SOS, Full-SOS (ms/degC)

edge_noise = @(im) max([reshape(im(1:10,:),[],1); reshape(im(end-9:end,:),[],1); reshape(im(:,1:10),[],1); reshape(im(:,end-9:end),[],1)]);
noise_sos = edge_noise(double(S1f.img4isos(:,:,PR_TP,1)));
[ri,ci] = ndgrid(1:100,1:100);
oil100 = logical(S3f.yangaib(sub2ind([128 128], round(1+(ri-1)*127/99), round(1+(ci-1)*127/99))));
oil100 = circshift(oil100,g.oil_from_legacy);
rod = g.rod_100;
nearrod = false(100); nearrod(max(rod(1)-3,1):min(rod(2)+3,100), max(rod(3)-3,1):min(rod(4)+3,100)) = true;
box_r = unique(round(1+((51:90)-1)*99/127)); box_c = unique(round(1+((45:50)-1)*99/127));
rng(7,'twister');

variants = {'PFA',1,g.pfa_roi; 'Full',3,g.full_roi};
OUT = struct();
for v = 1:2
    vn = variants{v,1}; vi = variants{v,2}; roi = variants{v,3};
    d = dir(fullfile(prdir,sprintf('%s_tp%d_chunk*.mat',vn,PR_TP)));
    d = d(~contains({d.name},'chunk99'));
    sos3 = []; mip3 = []; T2s = []; T2m = [];
    for j = 1:numel(d)
        X = load(fullfile(prdir,d(j).name));
        sos3 = cat(4,sos3,double(X.sos3)); mip3 = cat(4,mip3,double(X.mip3));
        T2s = cat(3,T2s,double(X.T2s)); T2m = cat(3,T2m,double(X.T2m));
    end
    % Use the same number of replicas for both variants (the first PR_NR_MAX).
    if size(T2s,3) < PR_NR_MAX
        error('rev:TooFewReplicas', '%s: %d replicas found, %d needed. Fix: run rev_pseudoreplica.m for more chunks (PR_CHUNK = 1, 2, ...).', vn, size(T2s,3), PR_NR_MAX);
    end
    sos3 = sos3(:,:,:,1:PR_NR_MAX); mip3 = mip3(:,:,:,1:PR_NR_MAX);
    T2s = T2s(:,:,1:PR_NR_MAX); T2m = T2m(:,:,1:PR_NR_MAX);
    NR = size(T2s,3);
    T2b = double(S1f.T2MAPsos(:,:,PR_TP,vi));
    mask = T2b>0 & T2b<=200 & double(S1f.img4isos(:,:,PR_TP,vi))>=3*noise_sos & ~oil100 & ~nearrod ...
        & all(T2s>0 & T2s<400,3) & all(T2m>0 & T2m<400,3);
    fprintf('\n===== %s: %d replicas, %d analysis voxels =====\n', vn, NR, nnz(mask));
    I16 = squeeze(sos3(:,:,1,:));
    maps = {I16, T2s, T2m}; mlab = {'SOS image TE16','T2 SOS','T2 MIP'};
    for i = 1:3
        [ro,pe,dg] = lagz(maps{i},mask);
        fprintf('%-15s lag corr RO 1-3: %+.3f %+.3f %+.3f | PE 1-5: %+.3f %+.3f %+.3f %+.3f %+.3f | diag %+.3f\n', mlab{i}, ro, pe, dg);
        OUT(v).lag(i).name = mlab{i}; OUT(v).lag(i).ro = ro; OUT(v).lag(i).pe = pe; OUT(v).lag(i).diag = dg;
    end
    idx = find(mask); cte = zeros(1,3); pairs = [1 2;1 3;2 3];
    for p = 1:3
        a = reshape(sos3(:,:,pairs(p,1),:),[],NR); b = reshape(sos3(:,:,pairs(p,2),:),[],NR);
        a = a(idx,:) - mean(a(idx,:),2); b = b(idx,:) - mean(b(idx,:),2);
        cte(p) = mean(sum(a.*b,2)./sqrt(sum(a.^2,2).*sum(b.^2,2)));
    end
    fprintf('cross-TE same-voxel corr (SOS): TE16-24 %+.3f TE16-32 %+.3f TE24-32 %+.3f\n', cte);
    OUT(v).cross_te = cte;
    a = I16(:,1:end-1,:); b = squeeze(sos3(:,2:end,3,:)); mm = mask(:,1:end-1) & mask(:,2:end);
    a = reshape(a,[],NR); b = reshape(b,[],NR); a = a(mm(:),:) - mean(a(mm(:),:),2); b = b(mm(:),:) - mean(b(mm(:),:),2);
    OUT(v).cross_te_pe1 = mean(sum(a.*b,2)./sqrt(sum(a.^2,2).*sum(b.^2,2)));
    fprintf('cross-TE PE-neighbour corr TE16(x) vs TE32(x+1 PE): %+.3f\n', OUT(v).cross_te_pe1);

    for c = 1:2
        TT = {T2s,T2m}; T = TT{c}; clab = {'SOS','MIP'};
        for k = 1:4
            rr = roi(k,1)-1:roi(k,1)+1; cc = roi(k,2)-1:roi(k,2)+1;
            X = reshape(T(rr,cc,:),9,NR)';
            [rbar,sdm,sdn] = roistats(X);
            OUT(v).roi(c,k) = struct('rbar',rbar,'sd_mean',sdm,'sd_naive',sdn, ...
                'sd_vox',sqrt(mean(var(X))),'baseline_T2',mean(reshape(T2b(rr,cc),[],1)));
        end
        Q = OUT(v).roi(c,:);
        fprintf('%s %s ROI I-IV: mean off-diag r %s | SD(ROI-mean T2) %s ms vs naive SD/3 %s ms (ratio %s) | voxel SD %s ms | baseline T2 %s ms\n', ...
            vn, clab{c}, mat2str([Q.rbar],2), mat2str([Q.sd_mean],3), ...
            mat2str([Q.sd_naive],3), mat2str([Q.sd_mean]./[Q.sd_naive],3), mat2str([Q.sd_vox],3), mat2str([Q.baseline_T2],3));
    end
    Q = OUT(v).roi(1,:);
    n_tp_ = size(S1f.T2MAPsos,3);
    fprintf('%s SOS noise SD of ROI-mean dT (deviation from %d-point mean, / pooled slope %.4f): %s degC\n', vn, n_tp_, m_pool(v), ...
        mat2str([Q.sd_mean]*sqrt((n_tp_-1)/n_tp_)/m_pool(v),3));

    sigS = sqrt(mean(reshape(var(I16(box_r,box_c,:),0,3),[],1)));
    M16 = squeeze(mip3(:,:,1,:)); sigM = sqrt(mean(reshape(var(M16(box_r,box_c,:),0,3),[],1)));
    snrS = mean(S1f.img4isos(box_r,box_c,PR_TP,vi),'all')/sigS; snrM = mean(S1f.img4i(box_r,box_c,PR_TP,vi),'all')/sigM;
    fprintf('%s replica SNR at time point %d (TE16, reference rectangle %d voxels): SOS %.1f  MIP %.1f\n', vn, PR_TP, numel(box_r)*numel(box_c), snrS, snrM);
    OUT(v).snr = [snrS snrM]; OUT(v).sigma = [sigS sigM]; OUT(v).NR = NR; OUT(v).nmask = nnz(mask); OUT(v).tp = PR_TP;
end

fprintf('\nLINEAR VIEW-SHARING MODEL (Full: central ky -5..4 TE-specific; PFA approximated without POCS: central ky -3..2)\n');
for T2 = [50 70 90]
    S = exp(-(16:8:32)/T2)'; gg = zeros(3,1); h = 1e-6;
    f0 = nlsT2(S);
    for n = 1:3, Sp = S; Sp(n) = Sp(n)+h; gg(n) = (nlsT2(Sp)-f0)/h; end
    ky = -50:49;
    for vv = 1:2
        if vv == 1, cen = -5:4; acq = ky; else, cen = -3:2; acq = -10:49; end
        P = zeros(size(ky)); P(ismember(ky,acq)) = sum(gg)^2; P(ismember(ky,cen)) = sum(gg.^2);
        rho = arrayfun(@(L) sum(P.*cos(2*pi*ky*L/100))/sum(P), 1:5);
        fprintf('T2 %d ms %-4s: sum g^2 / (sum g)^2 = %.1f; predicted PE corr lag 1-5: %s\n', T2, tern(vv==1,'Full','PFA*'), ...
            sum(gg.^2)/sum(gg)^2, mat2str(rho,3));
    end
end
save(fullfile(outdir,'pseudoreplica_summary.mat'),'OUT');

function [ro,pe,dg] = lagz(M,mask)
Z = (M - mean(M,3))./std(M,0,3);
ro = arrayfun(@(L) pz(Z,mask,L,0),1:3); pe = arrayfun(@(L) pz(Z,mask,0,L),1:5); dg = pz(Z,mask,1,1);
end
function r = pz(Z,mask,dr,dc)
[n1,n2,~] = size(Z);
m = mask(1:n1-dr,1:n2-dc) & mask(1+dr:n1,1+dc:n2);
a = Z(1:n1-dr,1:n2-dc,:); b = Z(1+dr:n1,1+dc:n2,:);
p = mean(a.*b,3)*size(Z,3)/(size(Z,3)-1);
r = mean(p(m));
end
function [rbar,sdm,sdn] = roistats(X)
S = cov(X); o = ones(9,1);
Rm = corrcov(S); rbar = mean(Rm(~eye(9)));
sdm = sqrt(o'*S*o)/9; sdn = sqrt(trace(S)/9)/3;
end
function T2 = nlsT2(S)
n = (1:3)';
J = @(t) sum(S.^2) - (S'*exp(-n*8/t))^2/sum(exp(-2*n*8/t));
T2 = fminbnd(J,5,400,optimset('TolX',1e-10));
end
function s = tern(c,a,b), if c, s=a; else, s=b; end, end
