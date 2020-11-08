% G = nuft_kern(kx, ky, res, w)
%
% Function that computes the convolution kernel related to the operation
% of the k-space measurement of a 2D image followed by its adjoint.
% 
%                    Ns
%       G[px,py] =   sum  w[n]*exp(-j*2*pi*(kx[n]*px/FOVx+ky[n]*py/FOVy)),
%                    n=1
% with -FOVx+1 <= px <= FOVx-1 and -FOVy+1 <= py <= FOVy-1.
%
% Inputs: the 2D k-space trajectory, the size of FOV and the optional weighting vector.
%
% Output: the 2D kernel image.
%
% No gridding is involved. Note that FFT-based implementations can be much faster.
% The implementation aims at minimizing the amount of trigonometric function calls.
% Observed x40 speedup compared to the trivial MATLAB implementation (code below).
%
% Copyright, Matthieu Guerquin-Kern, 2012

function G = nuft_kern(kx, ky, res, w)

G = zeros(2*res-1);
[X,Y] = ndgrid((-res(1)+1:res(1)-1)/res(1),(-res(2)+1:res(2)-1)/res(2));
for i=1:length(kx)
	G = G+w(i)*exp(-1i*2*pi*(kx(i)*X+ky(i)*Y));
end
