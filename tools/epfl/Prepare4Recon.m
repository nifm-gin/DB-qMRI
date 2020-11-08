%% [a,A,P] = Prepare4Recon(m, k, sens, support, V)
%
% Prepares data for reconstruction
%
% Copyright, Matthieu Guerquin-Kern, 2012

function [a,A,P] = Prepare4Recon(m, k, sens, support, V)

%% Initializations
if nargin==5
	param.V = V;
end
if nargin>3
	param.sensitivity = sens.*repmat(support,[1,1,size(sens,3)]);
else
	param.sensitivity = sens;
    support = ones(size(sens(:,:,1)));
end

%% Computations
rsos = sqrt(sum(abs(sens).^2,3));
param.k = k;
param.res = size(sens(:,:,1));
param.method = 'nuft gg fast';

a = EH(m,param); % backproject the data
[y, param] = EHVE(a, param); % this way we precompute the kernel inside param
A = @(x) EHVE(x,param); % create the function handle
P = ones(param.res);
P(support) = 1./(rsos(support)+1e-5*max(rsos(support(:))));
