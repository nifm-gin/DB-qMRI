function m = nuft_gg_forw(x,st)
% This routine performs forward non-uniform Fourier transform using FFT
% and a Gaussian interpolation kernel. It implements Greengard's fast
% Gaussian gridding.
%
% The input structure must be initialized with nuft_gg_init.m
%
% Copyright, Matthieu Guerquin-Kern, 2012

% Deconvolve, Zero-pad and compute DFT
C =  zeros(st.Mr);
C(1:st.M(1),1:st.M(2)) = st.E4.*x;
C = circshift(C,-floor(st.M/2));
C = fftn(C)/prod(st.Mr);

% Perform fast gridding and postfilter
m = st.E1.*nuft_forw_gridding(C,st.m1,st.m2,st.E2x,st.E2y,st.E3x,st.E3y);