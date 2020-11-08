%% [x, t, delta] = ReconTV(a, A, lambda, x0, Nit, Ninter, P)
%
% Performs Total Variation penalized reconstruction trying to solve the 
% linear system
%                       Ax = a,
% with A a Hermitian symmetric and positive-definite matrix.
%
% Optionally, an image reprensenting a diagonal matrix P can be used for
% preconditionning such that the following system is solved
%                       Mu = b,
% with M = PAP, x = Pu and b = Pa.
%
% The Total variation considered is isotropic and the algorithm used to
% perform the reconstruction is Iteratively Reweighted Least-Squares [1].
%
% [1] B. Wohlberg and P. Rodríguez, "An iteratively reweighted norm
% algorithm for minimization of total variation functionals,"
% IEEE Signal Processing Letters, vol. 14, no. 12, pp. 948-951, Dec. 2007.
%
% Copyright, Matthieu Guerquin-Kern, 2012

function [x,t,delta] = ReconTV(a, A, lambda, varargin)

%% Default inputs
if nargin<3
    error('ReconTV:TooFewInputs', ...
        'requires at least 3 inputs');
end
numvarargs = length(varargin);
if numvarargs > 4
    error('ReconTV:TooManyInputs', ...
        'requires at most 7 inputs');
end
optargs = {zeros(size(a)) 30 10 ones(size(a))};
optargs(1:numvarargs) = varargin;
[x0,  Nit, Ninter, P] = optargs{:};
t = 0;
delta = [];
h = waitbar(0,'TV-IRLS reconstruction...');

%% Computations
finitediff = @(x,dir) x-circshift(x,dir);
DTwD = @(x,w) finitediff(w.*finitediff(x,[1 0]),[-1 0]) + finitediff(w.*finitediff(x,[0 1]),[0 -1]);
x = x0;
for i=1:Nit
	w = (abs(finitediff(x,[1 0])).^2+abs(finitediff(x,[0 1])).^2);
    w = .5*lambda./sqrt(w + 1e-12*max(w(:)));
	M = @(x) A(x)+DTwD(x,w);
	[x,tmp,dmp] = ReconCG(a,M,x,Ninter,P);
    t = [t; t(end)+tmp(2:end)];
    delta = [delta; dmp(2:end)];
    waitbar(i/Nit,h);
end
t = t(2:end);
close(h);
