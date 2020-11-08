%% V = EstimateCovarianceMatrix(n)
%
% Estimates the cross-channel covariance matrix out of noise only data.
% Returns the More-Penrose pseudoinverse of the covariance matrix.
%
% INPUTS:    * a [N*Nc] complex-valued array
%
% OUTPUTS:   * a [Nc*Nc] Hermitian symmetric array
%
% Copyright, Matthieu Guerquin-Kern, 2012

function V = EstimateCovarianceMatrix(n)

n = n-ones(size(n,1),1)*mean(n,1);
V = pinv(n'*n/size(n,1));
