function [xs_hat, xw_hat] = WavAnalysisSep(xs_hat, htldstr_hat, gtldstr_hat, decimation)
% [xs_hat, xw_hat] = WavAnalysisSep(x_hat, htldstr_hat, gtldstr_hat, decimation)
%
% Frequency-domain implementation of the analysis side of a
% (dyadically decimated or undecimated) separable filter bank.
% Generally, this is a 2^d-channel filter bank, where d is the
% dimensionality of the input signal (see also the note below).
%
% x_hat:       DFT of the signal to be decomposed.
% htldstr_hat,
% gtldstr_hat: Complex conjugates of the DFTs of the (scaling-function
%              and wavelet) analysis filters (cell arrays corresponding
%              to the dimensionality of x_hat).
% decimation:  True or false.
%
% NOTE: if the cell of htldstr corresponding to the d-th dimension is
%       empty, the analysis will not be performed along that dimension.
%       This can be used to obtain dimension-dependent decomposition depths.
%
% xs_hat:      DFT of the coarsest-scale scaling function coefficients.
% xw_hat:      DFTs of the wavelet coefficients (cell array corresponding
%              to the wavelet subbands).
%
% (c) Cedric Vonesch, 2006.11.06-2007.10.21

% The dimensions to be decomposed
DecomposeDim = ~cellfun(@isempty, htldstr_hat);

% Allocate DFTs of the wavelet coefficients
xw_hat = cell(1, 2^sum(DecomposeDim)-1);

if length(decimation)==1
    decimation = [decimation,decimation]; % first counts for spline coeff, second for wavelet coeff
end

% Filter-bank analysis
ds = 0;
for d = 1:numel(DecomposeDim)
	if DecomposeDim(d)
		if nargout > 1
			xw_hat{2^ds} = WavAnalysis1D(xs_hat, gtldstr_hat{d}, d, decimation(2));
			for s = 1:2^ds-1
				xw_hat{2^ds+s} = WavAnalysis1D(xw_hat{s}, gtldstr_hat{d}, d, decimation(2));
				xw_hat{s} = WavAnalysis1D(xw_hat{s}, htldstr_hat{d}, d, decimation(2));
			end
		end
		xs_hat = WavAnalysis1D(xs_hat, htldstr_hat{d}, d, decimation(1));
		ds = ds + 1;
	end
end