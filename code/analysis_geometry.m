function g = analysis_geometry(R, reconstruction_frame)
%ANALYSIS_GEOMETRY One authority for native coordinates after the POCS fix.
% The archived MAT records refer to the legacy PFA image origin. Correcting
% the inverse Fourier transform moves PFA image content +1 native column.
% Map its ROI indices by the same amount to preserve the original samples.
% Full reconstruction and its actual statistical ROI indices are unchanged.
% No ROI is reselected and no mapping is chosen using temperature errors.
assert(strcmp(reconstruction_frame,'centered-fft-v1'), ...
    'vsmsi:LegacyFrame','Reconstruct maps with the corrected POCS before analysis.');
g.pfa_from_legacy = [0 1];
g.pfa_roi = double(R.ROI_2DVSMSI_100) + g.pfa_from_legacy;
g.full_roi = double(R.ROI_2DVSMSI_100);
g.prf_roi = double(R.ROI_PRF_128);
g.pfa_series = R.ROI_2DVSMSI_100_series;
g.full_series = R.ROI_2DVSMSI_100_series;
for k=1:numel(g.pfa_series)
    g.pfa_series{k} = double(g.pfa_series{k}) + g.pfa_from_legacy;
end
% Native pixel size, mm, [rows = readout, columns = phase encoding].
% The ADC samples 0.52 ms of the 0.54 ms readout plateau, so the readout FOV is
% 256*0.54/0.52 = 265.85 mm (calculateKspacePP on the scanned .seq files and
% oil-bottle positions against GRE).
% The phase-encoding FOV is the nominal 256 mm. GRE is 256 mm in both.
g.px_100 = [256*0.54/0.52, 256]/100;
g.px_128 = [2 2];
% Isocentre index of the centred-FFT grids (N/2+1).
g.centre_100 = 51;
g.centre_128 = 65;
% [first row, last row, first column, last column]
g.rod_100 = [40 70 49 52];
g.rod_128 = [51 91 62 67];
% Integer display approximation of the native 100-grid rod after resampling
% to the GRE grid with the physical mapping above.
g.rod_resampled_128 = round(g.centre_128 + ...
    (g.rod_100 - g.centre_100).*g.px_100([1 1 2 2])./g.px_128([1 1 2 2]));
g.oil_from_legacy = [0 1];
g.display_size = 128;
end
