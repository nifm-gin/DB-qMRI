%% [y, param] = EHVE(x, param)
%
% Function that performs the matrix-vector multiplication EH*V*E with E the
% matrix that models a 2D parallel MRI experiment for a given arbitrary k-space
% trajectory, E0H it Hermitian transpose, and V a noise-decorrelating matrix.
% This implementation is equivalent but faster than the consecutive use of
% functions E and EH.
%
% INPUTS:       * the original 2D image
%               * a structure containing the parameters (sensistivity maps, 
%			kernel, inverse covariance matrix,... ).
%
% OUTPUTS:      * the resulting 2D image
%               * (optionally) a structure carrying precomputed data to fasten 
%			future calls.
%
% SEE: E.m EH.m EHE0.m
%
% Copyright, Matthieu Guerquin-Kern, 2012

function [y, param] = EHVE(x, param)

%% Initializations
if ~isfield(param,'sensitivity')
    [y,param] = EHE0(x, param);
    return;
end
res = size(x);
y = zeros(res);
Nc = size(param.sensitivity,3);

%% Computations
if isfield(param,'V')
    data = zeros(prod(res),Nc);
    for ind_coil = 1:Nc
        [tmp, param] = EHE0(x.*param.sensitivity(:,:,ind_coil),param);
        data(:,ind_coil) = reshape(conj(param.sensitivity(:,:,ind_coil)).*tmp,[prod(res),1]);
    end
    y = reshape(sum(data*param.V,2),res);    
else
    for ind_coil = 1:Nc
        [tmp, param] = EHE0(x.*param.sensitivity(:,:,ind_coil),param);
        y = y + conj(param.sensitivity(:,:,ind_coil)).*tmp;
    end
end
