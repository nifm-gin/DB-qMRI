function [N, D, M] = GetDimensions(x)
% [N, D, M] = GetDimensions(x)
%
% Returns the size and dimension of a signal x.
%
% x:    The signal.
%
% N:    The dimensions of x.
% D:    The number of dimensions of x (length of N).
% M:    Corresponding number of wavelet subbands per
%       scale for a separable wavelet transform.
%
% Note: this function correctly returns D = 1 if x
% is a 1D signal, unlike the Matlab routine ndims().
%
% (c) Cedric Vonesch, 2007.03.17-2008.04.03

N = size(x);
D = length(N);
if D == 2 && any(N == 1)
	D = 1;
	N = prod(N);
end
M = 2^D-1;