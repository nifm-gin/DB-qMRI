function [h_hat, g_hat, htld_hat, gtld_hat] = BiorthWavFilters1D(N, family, j)
% [h_hat, g_hat, htld_hat, gtld_hat] = BiorthWavFilters1D(N, family, j)
%
% Generate the DFTs of the scaling and wavelet filters of
% various orthonormal or biorthogonal 1D wavelet families.
% In the orthonormal case:
% - the filters are normalized to unit l2 norm;
% - the analysis and synthesis filters are the same.
%
% N:        Filter length.
% family:   Wavelet family (see code).
% j:        Optional scale argument (for non-stationary
%           wavelet bases).
%
% h_hat:    DFT of the synthesis scaling filter.
% g_hat:    DFT of the synthesis wavelet filter.
% htld_hat: DFT of the analysis scaling filter.
% gtld_hat: DFT of the analysis wavelet filter.
%
% (c) Cedric Vonesch, 2006.09.05-2008.04.08

% Scaling filter
switch lower(family)
	case {'haar', 'daub1'} % Haar
		h = zeros(N, 1);
		h([end, 1]) = [1, 1]/sqrt(2);
		h_hat = fft(h);
	case 'daub2' % 2nd-order Daubechies minimal phase
		h = zeros(N, 1);
		h([end-1, end, 1, 2]) = [1+sqrt(3), 3+sqrt(3), 3-sqrt(3), 1-sqrt(3)]/4/sqrt(2);
		h_hat = fft(h);
    case 'daub4'
        h = zeros(N, 1);
        h([end-3:end, 1:4]) = [0.230377813308896   0.714846570552915   0.630880767929859  -0.027983769416859  -0.187034811719093   0.030841381835561   0.032883011666885  -0.010597401785069];
        h_hat = fft(h);
	case 'sym4' % 4th-order symlet
		h = zeros(N, 1);
		h([end-3:end, 1:4]) = wfilters('sym4');
		h_hat = fft(h);
	case {'spline1','spline2','spline3','spline4','spline5','spline6'} % pth-order orthonormalized spline
		p = str2double(family(7));
		h_hat = (1 + exp(1i*2*pi*(0:N-1).'/N)).^p; % p-th order zero at pi
		h_hat = exp(1i*2*pi*(-p/2)*(0:N-1).'/N) .* h_hat; % Center the filter at 0
		h_hat = sqrt(2) * h_hat ./ sqrt(abs(h_hat).^2 + fftshift(abs(h_hat).^2)); % Orthonormalize the filter
	case 'sym8' % 8th-order symlet
		h = zeros(N, 1);
		h([end-7:end, 1:8]) = wfilters('sym8');
		h_hat = fft(h);
	case {'sinc', 'shannon'} % Sinc
		h_hat = [sqrt(2)*ones(1, ceil(N/4)), zeros(1, N/2), sqrt(2)*ones(1, floor(N/4))].';
	case {'97', '9/7'} % Biorthogonal 9/7-tap Daubechies
		h = zeros(N, 1);
		g = zeros(N, 1);
		htld = zeros(N, 1);
		gtld = zeros(N, 1);
		[htld([end-4:end, 1:5]), gtld([end-2:end, 1:7]), h([end-3:end, 1:6]), g([end-3:end, 1:6])] = wfilters('bior4.4');
		h_hat = fft(h);
		g_hat = fft(g);
		htld_hat = fft(htld);
		gtld_hat = fft(gtld);
		return;
	case 'espline' % Exponential B-spline, Ildar's code
		alpha0 = [-1.0204 -2.9697 -0.3247+0.5487*1i  -0.3247-0.5487*1i];
		gamma0 = [2.5524];
		[h_hat, g_hat, htld_hat, gtld_hat] = precompute_filters(N, 2^(j-1)*alpha0, 2^(j-1)*gamma0, 1);
		h_hat = h_hat.';
		g_hat = g_hat.';
		htld_hat = htld_hat.';
		gtld_hat = gtld_hat.';
		return;
	otherwise
		error('Unknown wavelet family.');
end

% Complete the filter set in the orthonormal case
g_hat = exp(-1i*2*pi/N*(0:N-1)).' .* conj(fftshift(h_hat));
htld_hat = h_hat;
gtld_hat = g_hat;