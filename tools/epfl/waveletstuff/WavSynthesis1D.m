function xs_hat = WavSynthesis1D(xs_hat, h_hat, d, decimation)
% x_hat = WavSynthesis1D(xs_hat, h_hat, d, decimation)
%
% Frequency-domain implementation of the synthesis side of a single
% channel within a (dyadically decimated or undecimated) filter bank,
% along the d-th dimension of a signal.
%
% xs_hat,
% xw_hat:     DFTs of the scaling function and wavelet coefficients.
% h_hat,
% g_hat:      DFTs of the (scaling function and wavelet) synthesis filters.
% d:          The dimension to be reconstructed.
% decimation: True or false.
%
% x_hat:      DFT of the reconstructed signal.
%
% (c) Cedric Vonesch, 2007.10.18-2008.04.04

if decimation
	xs_hat = UpSampleDim(xs_hat, d);
end
xs_hat = FilterSep(xs_hat, d, h_hat);