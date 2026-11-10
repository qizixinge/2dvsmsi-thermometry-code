%% Synthetic regression tests for pocsnew
% No study data are read. The tests exercise both missing sides and both
% spatial PF dimensions on even, non-square matrices.

testDir = fileparts(mfilename('fullpath'));
addpath(fullfile(fileparts(testDir),'code'));

cases = { ...
    struct('name','dim2-lower-paper-size', 'sz',[100 100], 'pfDim',2, 'Nsmp',60, 'isLower',true), ...
    struct('name','dim2-upper-paper-size', 'sz',[100 100], 'pfDim',2, 'Nsmp',60, 'isLower',false), ...
    struct('name','dim1-lower-nonsquare',  'sz',[96 80],   'pfDim',1, 'Nsmp',58, 'isLower',true), ...
    struct('name','dim1-upper-nonsquare',  'sz',[96 80],   'pfDim',1, 'Nsmp',58, 'isLower',false) ...
};

for c = 1:numel(cases)
    tc = cases{c};
    [imRef, kspRef] = make_phantom(tc.sz, tc.pfDim);
    kspPF = mask_pf(kspRef, tc.pfDim, tc.Nsmp, tc.isLower);

    [imPOCS, kspPOCS] = pocsnew(kspPF, 30, false);
    imZF = centred_ifft2(kspPF);

    assert(isequal(size(imPOCS), size(imRef)), '%s: image size changed', tc.name);
    assert(isequal(size(kspPOCS), size(kspRef)), '%s: k-space size changed', tc.name);
    assert(all(isfinite(imPOCS(:))) && all(isfinite(kspPOCS(:))), '%s: non-finite output', tc.name);

    % The returned image must be the standard centred inverse FFT of the
    % returned k-space. This detects the legacy one-pixel shift/phase ramp.
    imFromKspace = centred_ifft2(kspPOCS);
    pairError = norm(imPOCS(:)-imFromKspace(:)) / max(norm(imFromKspace(:)), eps);
    assert(pairError < 5e-12, '%s: im and kspFull are not a centred Fourier pair (%.3g)', ...
        tc.name, pairError);

    rmsePOCS = norm(abs(imPOCS(:))-abs(imRef(:))) / norm(abs(imRef(:)));
    rmseZF = norm(abs(imZF(:))-abs(imRef(:))) / norm(abs(imRef(:)));
    % Absolute bound AND a margin over the zero-filled baseline on the same data: plain zero-filling already
    % reaches 0.013-0.018 here, so a threshold above that would pass a broken or missing POCS step.
    assert(rmsePOCS < 0.015, '%s: POCS error %.4g above the absolute bound 0.015', tc.name, rmsePOCS);
    assert(rmsePOCS < 0.7*rmseZF, '%s: POCS error %.4g is not clearly below the zero-filled error %.4g', ...
        tc.name, rmsePOCS, rmseZF);

    % pocsnew deliberately blends Ntrans acquired lines at the PF boundary.
    % Test exact data consistency only in the documented strict region.
    centreLine = floor(tc.sz(tc.pfDim)/2) + 1;
    if tc.isLower && mod(tc.sz(tc.pfDim),2) == 0
        centreLine = centreLine - 1;
    end
    szSym = 2*(tc.Nsmp-centreLine)+1;
    Ntrans = floor((szSym-1)/3);
    if tc.isLower
        pure = (tc.sz(tc.pfDim)-tc.Nsmp+Ntrans+1):tc.sz(tc.pfDim);
    else
        pure = 1:(tc.Nsmp-Ntrans);
    end
    d = select_dim(kspPOCS-kspPF, tc.pfDim, pure);
    ref = select_dim(kspPF, tc.pfDim, pure);
    relConsistency = norm(d(:)) / max(norm(ref(:)), eps);
    assert(relConsistency < 5e-12, '%s: strict acquired samples changed (relative %.3g)', ...
        tc.name, relConsistency);

    fprintf('PASS %-27s POCS %.5f, zero-fill %.5f, pair %.3g, consistency %.3g\n', ...
        tc.name, rmsePOCS, rmseZF, pairError, relConsistency);
end

fprintf('PASS: %d reconstruction-level POCS cases.\n', numel(cases));

function [im, ksp] = make_phantom(sz, rampDim)
    [y,x] = ndgrid(linspace(-1,1,sz(1)), linspace(-1,1,sz(2)));
    mag = 0.9*exp(-4*(x.^2+1.3*y.^2)) ...
        + 0.45*exp(-35*((x+0.35).^2+(y-0.20).^2)) ...
        + 0.25*(abs(x-0.25)<0.16 & abs(y+0.30)<0.11);
    smoothPhase = 0.18*sin(pi*x).*cos(pi*y);
    if rampDim == 1
        smoothPhase = smoothPhase + 2*pi*2*(0:sz(1)-1)'/sz(1);
    else
        smoothPhase = smoothPhase + 2*pi*2*(0:sz(2)-1)/sz(2);
    end
    im = reshape(mag.*exp(1i*smoothPhase), [1 sz]);
    ksp = fftshift(fftshift(fft(fft(ifftshift(ifftshift(im,2),3),[],2),[],3),2),3);
end

function kspPF = mask_pf(ksp, pfDim, Nsmp, isLower)
    kspPF = ksp;
    N = size(ksp, pfDim+1);
    subs = {':',':',':'};
    if isLower
        subs{pfDim+1} = 1:(N-Nsmp);
    else
        subs{pfDim+1} = (Nsmp+1):N;
    end
    kspPF(subs{:}) = 0;
end

function im = centred_ifft2(ksp)
    im = fftshift(fftshift(ifft(ifft(ifftshift(ifftshift(ksp,2),3),[],2),[],3),2),3);
end

function out = select_dim(in, pfDim, idx)
    subs = {':',':',':'};
    subs{pfDim+1} = idx;
    out = in(subs{:});
end
