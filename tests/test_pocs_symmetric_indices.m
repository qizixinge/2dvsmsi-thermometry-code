%% Regression tests for the symmetric low-resolution POCS window
% No study data are read. Part 1 isolates the index convention used after pocsnew flips a
% lower-side acquisition so that acquired samples occupy 1:Nsmp (index arithmetic only).
% Part 2 calls pocsnew itself on synthetic lower-side partial-Fourier input and checks the
% flip and symmetric-window behaviour through its output.

testDir = fileparts(mfilename('fullpath'));
addpath(fullfile(fileparts(testDir),'code'));

cases = { ...
    struct('name','paper 100x100 lower side', 'sz',[100 100], 'pfDim',2, 'Nsmp',60, 'isLower',true,  'expectedCentre',50, 'expectedIdx',40:60), ...
    struct('name','100x100 upper side',       'sz',[100 100], 'pfDim',2, 'Nsmp',60, 'isLower',false, 'expectedCentre',51, 'expectedIdx',42:60), ...
    struct('name','odd length lower side',    'sz',[80 101],  'pfDim',2, 'Nsmp',61, 'isLower',true,  'expectedCentre',51, 'expectedIdx',41:61), ...
    struct('name','odd length upper side',    'sz',[80 101],  'pfDim',2, 'Nsmp',61, 'isLower',false, 'expectedCentre',51, 'expectedIdx',41:61), ...
    struct('name','PF along first spatial dimension', 'sz',[96 100], 'pfDim',1, 'Nsmp',58, 'isLower',true, 'expectedCentre',48, 'expectedIdx',38:58) ...
};

for i = 1:numel(cases)
    tc = cases{i};
    [centreLine, idxSym] = candidate_indices(tc.sz, tc.pfDim, tc.Nsmp, tc.isLower);
    assert(centreLine == tc.expectedCentre, '%s: wrong centre', tc.name);
    assert(isequal(idxSym, tc.expectedIdx), '%s: wrong symmetric interval', tc.name);
    assert(all(idxSym >= 1 & idxSym <= tc.Nsmp), '%s: interval includes an unacquired line', tc.name);
    assert(numel(idxSym) == 2 * (tc.Nsmp - centreLine) + 1, '%s: wrong interval length', tc.name);
    assert((centreLine - idxSym(1)) == (idxSym(end) - centreLine), '%s: interval is not symmetric', tc.name);
end

% Reproduce the paper acquisition from dataSOS/dataprocessing:
% Nyo=100, pff=0.6 gives physical ky=-10:49 and native indices 41:100.
Nyo = 100;
Ny = 60;
physicalKy = -floor(Ny - Nyo/2):(Nyo/2 - 1);
nativeIndex = physicalKy + Nyo/2 + 1;
assert(isequal(sort(nativeIndex), 41:100));

% pocsnew detects a lower-side acquisition and flips it. Flipped index j maps
% to original native index 101-j and therefore to ky=50-j.
Nsmp = numel(nativeIndex);
legacyCentre = Nyo/2 + 1;
legacyIdx = (2 * legacyCentre - Nsmp):Nsmp;
[candidateCentre, candidateIdx] = candidate_indices([100 Nyo], 2, Nsmp, true);
legacyKy = 50 - legacyIdx;
candidateKy = 50 - candidateIdx;
assert(isequal(legacyIdx, 42:60));
assert(isequal(legacyKy, 8:-1:-10));
assert(candidateCentre == 50);
assert(isequal(candidateIdx, 40:60));
assert(isequal(candidateKy, 10:-1:-10));

% For odd N, the hard-coded centre is a half-integer. The colon endpoints
% happen to be integers here, but the resulting even-length interval is
% centred between samples rather than on the unique DC sample.
legacyOddCentre = 101/2 + 1;
legacyOddIdx = (2 * legacyOddCentre - 61):61;
assert(legacyOddCentre == 51.5);
assert(isequal(legacyOddIdx, 42:61));
assert(mean([legacyOddIdx(1), legacyOddIdx(end)]) == legacyOddCentre);
assert(~ismember(41, legacyOddIdx));

fprintf('PASS: %d generic cases plus the paper acquisition convention.\n', numel(cases));

%% Part 2: pocsnew on synthetic lower-side partial-Fourier data
% (a) Flip equivariance, odd length (index flip == ky -> -ky exactly). pocsnew flips a lower-side
%     acquisition, works on it as an upper-side one and flips the result back, so the lower-side
%     result must equal the index-flipped result of the mirrored upper-side input.
szOdd = [80 101];  NsmpOdd = 61;
[~, kspRefOdd] = make_phantom(szOdd);
kL = mask_pf(kspRefOdd, NsmpOdd, true);
kU = flip(kL, 3);
assert(all(all(kU(1,1,NsmpOdd+1:end) == 0)) && any(any(kU(1,1,1:NsmpOdd) ~= 0)), 'mirrored input is not an upper-side acquisition');
[~, kOutL] = pocsnew(kL, 30, false);
[~, kOutU] = pocsnew(kU, 30, false);
flipErr = norm(kOutL(:) - reshape(flip(kOutU,3),[],1)) / norm(kOutL(:));
assert(flipErr < 1e-10, 'lower-side pocsnew output is not the flipped upper-side output (relative %.3g)', flipErr);

% (b) Paper acquisition (100 lines, 60 acquired at the lower side, ky=-10..+49): the phase of the
%     reconstruction must follow the object phase, which needs the symmetric low-resolution window to
%     sit on the acquired central lines. The error must also be clearly below zero-filling.
szPaper = [100 100];  NsmpPaper = 60;
[imRefP, kspRefP] = make_phantom(szPaper);
kPaper = mask_pf(kspRefP, NsmpPaper, true);
imP = pocsnew(kPaper, 30, false);
imZ = centred_ifft2(kPaper);
objMask = abs(imRefP) > 0.3*max(abs(imRefP(:)));
phaseErr = @(im) mean(abs(angle(im(objMask).*conj(imRefP(objMask)))));
phP = phaseErr(imP);  phZ = phaseErr(imZ);
errP = norm(abs(imP(:))-abs(imRefP(:)))/norm(abs(imRefP(:)));
errZ = norm(abs(imZ(:))-abs(imRefP(:)))/norm(abs(imRefP(:)));
assert(phP < 0.02 && phP < 0.5*phZ, 'pocsnew phase error %.4g rad (zero-fill %.4g rad) too large', phP, phZ);
assert(errP < 0.7*errZ, 'pocsnew magnitude error %.4g not clearly below zero-filled %.4g', errP, errZ);
fprintf('PASS: pocsnew lower-side flip equivariance (%.2g) and paper-size phase %.4g rad (zero-fill %.4g rad).\n', flipErr, phP, phZ);

function [centreLine, idxSym] = candidate_indices(sz, pfDim, Nsmp, isLower)
    N = sz(pfDim);
    centreLine = floor(N/2) + 1;
    if isLower && mod(N,2) == 0
        % pocsnew has already reversed this dimension. For an even matrix,
        % the fftshift DC sample moves from N/2+1 to N/2 after reversal.
        centreLine = centreLine - 1;
    end
    idxSym = (2 * centreLine - Nsmp):Nsmp;
end

function [im, ksp] = make_phantom(sz)
    [y,x] = ndgrid(linspace(-1,1,sz(1)), linspace(-1,1,sz(2)));
    mag = 0.9*exp(-4*(x.^2+1.3*y.^2)) ...
        + 0.45*exp(-35*((x+0.35).^2+(y-0.20).^2)) ...
        + 0.25*(abs(x-0.25)<0.16 & abs(y+0.30)<0.11);
    phase = 0.18*sin(pi*x).*cos(pi*y) + 2*pi*2*(0:sz(2)-1)/sz(2);
    im = reshape(mag.*exp(1i*phase), [1 sz]);
    ksp = fftshift(fftshift(fft(fft(ifftshift(ifftshift(im,2),3),[],2),[],3),2),3);
end

function kspPF = mask_pf(ksp, Nsmp, isLower)
    kspPF = ksp;  N = size(ksp,3);
    if isLower, kspPF(:,:,1:(N-Nsmp)) = 0; else, kspPF(:,:,(Nsmp+1):N) = 0; end
end

function im = centred_ifft2(ksp)
    im = fftshift(fftshift(ifft(ifft(ifftshift(ifftshift(ksp,2),3),[],2),[],3),2),3);
end
