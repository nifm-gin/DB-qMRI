% m = nuft_forw(I, kx, ky)
%
% Function that performs the k-space measurements of a 2D image
%
%                 FOVx-1  FOVy-1
%       m[n] =     sum     sum   I[px,py]*exp(j*2*pi*(kx[n]*px/FOVx+ky[n]*py/FOVy)),
%                 px=0    py=0
% with 0 <= n <= Ns-1.
%
% Inputs: the 2D image I, the 2D k-space trajectory: kx, ky
%
% Output: the vector of k-space measurements.
%
% This is a MEX-file for MATLAB.  
%
% No gridding is involved. Note that FFT-based implementations can be much faster.
% The implementation aims at minimizing the amount of trigonometric function calls.
% Observed x20 speedup compared to the trivial MATLAB implementation (code below).
%
% Copyright, Matthieu Guerquin-Kern, 2012

function m = nuft_forw(I, kx, ky)

m = zeros(1,length(kx));
[X,Y] = ndgrid((0:size(I,1)-1)/size(I,1),(0:size(I,2)-1)/size(I,2));
for i=1:length(kx)
	m(i) = I(:).'*exp(1i*2*pi*(kx(i)*X(:)+ky(i)*Y(:)));
end
