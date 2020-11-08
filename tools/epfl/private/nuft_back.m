% I = nuft_back(m, kx, ky)
%
% Function that performs the adjoint operation of the k-space measurement of a 2D image
%
%                    Ns
%       I[px,py] =   sum  m[n]*exp(j*2*pi*(kx[n]*px/FOVx+ky[n]*py/FOVy)),
%                    n=1
% with 0 <= px <= FOVx-1, 0 <= py <= FOVy-1.
%
% Inputs: the k-space measurements m, the 2D k-space trajectory: kx,ky,
%		and the size of the output image FOVx,FOVy.
%
% Output: a 2D image
%
% This is a MEX-file for MATLAB.  
%
% No gridding is involved. Note that FFT-based implementations can be much faster.
% The implementation aims at minimizing the amount of trigonometric function calls.
% Observed x20 speedup compared to the trivial MATLAB implementation (code below).
%
% Copyright, Matthieu Guerquin-Kern, 2012

function I = nuft_back(m, kx, ky, s)

I = zeros(s);
[X,Y] = ndgrid((0:s(1)-1)/s(1),(0:s(2)-1)/s(2));
for i=1:length(kx)
	I = I+m(i).'*exp(-1i*2*pi*(kx(i)*X+ky(i)*Y));
end
