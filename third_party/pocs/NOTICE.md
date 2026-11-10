# POCS: original and modified version

`pocs.m` in this folder is the original, unmodified implementation of
Michael Völker, "MRI Partial Fourier reconstruction with POCS", MATLAB File
Exchange (https://www.mathworks.com/matlabcentral/fileexchange/39350-mri-partial-fourier-reconstruction-with-pocs),
Copyright (c) 2012, Michael Völker, distributed under the BSD 2-clause licence
in `LICENSE.txt` (this folder). It is kept for reference and is not called by
the analysis.

`code/pocsnew.m` is a modified version of that file, used by the analysis.
It is a derivative work and remains under the same BSD licence (the copyright
notice and conditions of `LICENSE.txt` apply to it). Modifications, in plain words:

1. The centre line of the symmetric acquired window is fixed at the geometric
   centre (DC position) of k-space in the partial-Fourier dimension, instead of
   being taken from the line of maximum signal. This matters for spectral bins
   with a non-zero frequency offset, where off-resonance moves the signal maximum.
2. The symmetric window is checked: the code stops with an error if a window
   symmetric about that centre cannot be formed.
3. The internal circular shift applied before the iterations is stored as an
   exact integer shift and is reversed exactly (also for odd matrix sizes), so
   the returned k-space is in the caller's original centred coordinates.
4. The image is obtained from the returned k-space with a centred inverse FFT,
   so that image and k-space are an exact Fourier pair for either missing side
   (the original flips an already transformed image).
5. Console progress messages are switched off.

The original code is credited in the header of `pocs.m` (Martin Blaimer, Uvo
Hoelscher, Michael Völker). The header of `code/pocsnew.m` carries the same credits.
