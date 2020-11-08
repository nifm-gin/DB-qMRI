%% [x, t, delta] = ReconCG(a, A, x0, Nit, P)
%
% Performs Conjugate gradient [1] reconstruction solving the linear system
%                       Ax = a,
% with A a Hermitian symmetric and positive-definite matrix.
%
% Optionally, a real-valued image representing a diagonal matrix P can be
% used for preconditionning such that the following system is solved
%                       Mu = b,
% with M = PAP, x = Pu and b = Pa.
%
% [1] K. P. Pruessmann, M. Weiger, P. Börnert, and P. Boesiger,
% "Advances in sensitivity encoding with arbitrary k-space trajectories,"
% Magnetic Resonance in Medicine, vol. 46, no. 4, pp. 638-651, 2001
%
% Copyright, Matthieu Guerquin-Kern, 2012

function [x,t,delta] = ReconCG(a, A, varargin)
%% Default inputs
if nargin<2
    error('ReconCG:TooFewInputs', ...
        'requires at least 2 inputs');
end
numvarargs = length(varargin);
if numvarargs > 3
    error('ReconCG:TooManyInputs', ...
        'requires at most 5 inputs');
end
optargs = {a 30 ones(size(a))};
optargs(1:numvarargs) = varargin;
[x0, Nit, P] = optargs{:};
t = zeros(Nit+1,1);
delta = zeros(Nit+1,1);
h = waitbar(0,'conjugate gradient reconstruction...');

%% Computations
Na2 = a(:)'*a(:);
a = P.*a;
u = x0./P;
r = a - P.*A(P.*u);
p = r;
delta(1) = r(:)'*r(:)/Na2;
for i=1:Nit
    t0 = clock();
	q = P.*A(P.*p);
	alpha = r(:)'*r(:)/(p(:)'*q(:));
	u = u + alpha*p;
	rnew = r - alpha*q;
	p = rnew + rnew(:)'*rnew(:)/(r(:)'*r(:))*p;
	r = rnew;
    t(i+1) = t(i) + etime(clock(),t0);
    delta(i+1) = r(:)'*r(:)/Na2;
    waitbar(i/Nit,h);
end
x = P.*u;
close(h);
