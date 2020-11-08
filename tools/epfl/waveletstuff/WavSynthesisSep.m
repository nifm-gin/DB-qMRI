function xs_hat = WavSynthesisSep(xs_hat, xw_hat, h_hat, g_hat, decimation)
% x_hat = WavSynthesisSep(xs_hat, xw_hat, h_hat, g_hat, decimation)
%
% Frequency-domain implementation of the synthesis side of a (dyadically
% decimated or undecimated) separable filter bank. Generally, this is a
% 2^d-channel filter bank, where d is the dimensionality of the input
% signal (see also the note in the documentation of WavAnalysisSep()).
%
% xs_hat:     DFT of the coarsest-scale scaling function coefficients.
% xw_hat:     DFTs of the wavelet coefficients (cell array corresponding
%             to the wavelet subbands).
% h_hat,
% g_hat:      DFTs of the (scaling-function and wavelet) synthesis filters,
%             (cell arrays corresponding to the dimensionality of x_hat).
% decimation: True or false.
%
% x_hat:      DFT of the reconstructed signal.
%
% (c) Cedric Vonesch, 2006.11.11-2008.04.04

% The dimensions do be reconstructed
ReconstructDim = ~cellfun(@isempty, h_hat);

if length(decimation)==1
    decimation = [decimation,decimation]; % first counts for spline coeff, second for wavelet coeff
end

% Filter-bank synthesis
ds = sum(ReconstructDim) - 1;
for d = numel(ReconstructDim):-1:1
	if ReconstructDim(d)
		if isempty(xs_hat)
			xs_hat = 0;
		else
			xs_hat = WavSynthesis1D(xs_hat, h_hat{d}, d, decimation(1));
		end
		if ~isempty(xw_hat)
			xs_hat = xs_hat + WavSynthesis1D(xw_hat{2^ds}, g_hat{d}, d, decimation(2));
			xw_hat{2^ds} = [];
			for s = 1:2^ds-1
				xw_hat{s} = WavSynthesis1D(xw_hat{s}, h_hat{d}, d, decimation(2)) + WavSynthesis1D(xw_hat{2^ds+s}, g_hat{d}, d, decimation(2));
				xw_hat{2^ds+s} = [];
			end
		end
		ds = ds - 1;
	end
end