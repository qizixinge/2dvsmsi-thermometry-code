function g = gre_params()
%GRE_PARAMS PRF-thermometry constants of the dual-echo GRE (Siemens gre_field_mapping) scan.
% DICOM header of the GRE series (e.g. series 0015/0016): ASCCONV alTE[0] = 5000 us,
% alTE[1] = 7460 us.  Series 0016 is the magnitude image (EchoNumber 1, TE 5 ms); the phase
% series 0015, 0017, ... are ImageType P with EchoNumber 2 / EchoTime 7.46 ms, i.e. the
% scanner's phase-DIFFERENCE image between echo 2 and echo 1.  The phase accrued between the
% two echoes therefore corresponds to dTE = 7.46 - 5.00 = 2.46 ms (not 5 ms).
%   g.B0     actual static field (T) = ImagingFrequency / gamma-bar. DICOM header (0018,0084) of the GRE series
%            0015/0016: ImagingFrequency 123.251054 MHz (the twix lFrequency, 123.251012 MHz, agrees to 0.3 ppm);
%            gamma-bar = 42.577478 MHz/T, so B0 = 2.8947 T, not the nominal 3 T (DICOM FieldStrength tag). The
%            scanner's own flMagneticFieldStrength entry (2.8936 T) is not used.
%   g.gamma  gyromagnetic ratio (rad/s/T)
%   g.dTE    effective echo-time difference (s)
%   g.window_s  data-acquisition duration (s): DICOM private header SliceMeasurementDuration = 25400 ms (TA 25.40 s)
%   g.scan_s    protocol total scan time (s): ASCCONV lTotalScanTimeSec = 27 (includes preparation cycles)
%   g.lam    B0*gamma*dTE (rad per unit chemical-shift change, i.e. phase = -lam * delta_s)
g.gamma_bar_MHz_T = 42.577478;                  % proton gyromagnetic ratio / 2pi, MHz/T
g.f0_MHz = 123.251054;                          % ImagingFrequency, DICOM tag (0018,0084)
g.B0 = g.f0_MHz/g.gamma_bar_MHz_T;
g.gamma = 2*pi*g.gamma_bar_MHz_T*1e6;
g.dTE = (7460 - 5000)*1e-6;
g.lam = g.B0*g.gamma*g.dTE;
g.window_s = 25.4;
g.scan_s = 27;
end
