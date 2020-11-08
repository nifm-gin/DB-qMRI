%% E.m
%
% Function that performs the matrix-vector multiplication that models a 2D
% parallel MRI experiment for given arbitrary k-space trajectory and
% sensitivity maps.
%
% INPUTS:       * the original 2D image
%               * a structure containing the parameters (k-space trajectory, 
%			sensitivity maps, method for computation).
%
% OUTPUTS:      * the k-space data vector
%               * (optionally) a structure carrying precomputed data to fasten future calls.
%
% SEE: EH.m EHVE.m
%
% Copyright, Matthieu Guerquin-Kern, 2012

function [m,param] = E(x, param)
%% Initializations

if ~isfield(param,'sensitivity')
    [m,param] = E0(x, param);
    return;
end

Nc = size(param.sensitivity,3);

%% Computation of the k-space measurements for each coil
for ind_coil = 1:Nc
    [tmp,param] = E0(x.*param.sensitivity(:,:,ind_coil), param);
    m(:,ind_coil) = tmp;
end
