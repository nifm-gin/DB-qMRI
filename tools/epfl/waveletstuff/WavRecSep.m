function xs_hat = WavRecSep(xs_hat, xw_hat, h_hat, g_hat, decimation)
% x_hat = WavRecSep(xs_hat, xw_hat, h_hat, g_hat, decimation)
%
% Frequency-domain implementation of a separable wavelet reconstruction.
% In particular:
% - the input must be the DFTs of the wavelet coefficients in each subband;
% - the function returns the DFT of the reconstructed signal.
%
% xs_hat:     DFT of the coarsest-scale scaling-function coefficients.
% xw_hat:     DFTs of the wavelet coefficients (two-level cell
%             array corresponding to every scale and every
%             wavelet subband).
% h_hat,
% g_hat:      DFTs of the (scaling function and wavelet)
%             synthesis filters (2D cell arrays corresponding
%             to every scale and every dimension).
% decimation: True or false.
%
% x_hat:  DFT of the reconstructed signal.
%
% (c) Cedric Vonesch, 2006.11.11-2008.04.04

% Decomposition depth
jmax = size(h_hat, 1);

% Reconstruction
for j = jmax:-1:1
	xs_hat = WavSynthesisSep(xs_hat, xw_hat{j}, h_hat(j, :), g_hat(j, :), decimation);
end