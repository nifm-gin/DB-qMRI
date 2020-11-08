%% [x,param] = EH(m, param)
%
% Function that performs the matrix-vector multiplication with the Hermitian
% transpose of the matrix that models a 2D MRI experiment for given 
% k-space trajectory and sensitivity maps.
%
% INPUTS:       * the k-space data vector
%               * a structure containing the parameters (k-space trajectory,
% 			method for computation, ...).
%
% OUTPUTS:      * the backprojected 2D image
%               * (optionally) a structure carrying precomputed data to fasten
% 			future calls.
%
% SEE: E.m EHVE.m
%
% Copyright, Matthieu Guerquin-Kern, 2012

function [x,param] = EH(m, param)
%% Initializations

if ~isfield(param,'sensitivity')
    [x,param] = EH0(m, param);
    return;
end

Nc = size(param.sensitivity,3);

%% Computation of the k-space measurements for each coil
x = zeros(param.res);
for ind_coil = 1:Nc
    [tmp,param] = EH0(m(:,ind_coil), param);
    x = x + conj(param.sensitivity(:,:,ind_coil)).*tmp;
end
