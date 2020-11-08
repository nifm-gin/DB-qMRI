function y = nuft_gg_back(m,st)
% This routine performs backward non-uniform Fourier transform using FFT
% and a Gaussian interpolation kernel. It implements Greengard's fast
% Gaussian gridding.
%
% The input structure must be initialized with nuft_gg_init.m
%
% Copyright, Matthieu Guerquin-Kern, 2012

% Prefilter the data and perform fast gridding
C = nuft_back_gridding(m.*conj(st.E1),st.m1,st.m2,st.Mr,st.E2x,st.E2y,st.E3x,st.E3y);

% Take inverse fft of C
c = ifftn(C);

% Get the relevant coefficients and deconvolve
c = circshift(c,floor(st.M/2));
y = st.E4.*c(1:st.M(1),1:st.M(2));