%% GenerateFullCartGrid.m
%
% Function that generates a 2D full Cartesian grid.
%
% INPUT:    * res : (2,1) vector or scalar that defines resolution
%
% OUTPUT:   * X1 : matrix where represents the first dimension
%           * X2 : matrix that represents the second definition
% 
% Note: N = nb of samples = prod(res)
%
% Matthieu Guerquin-Kern, Biomedical Imaging Group / EPF Lausanne,
% 30-10-2009 (dd-mm-yyyy)

function [X1,X2] = GenerateFullCart2DGrid(res)

if numel(res)==1
    res = res*[1,1];
end

[X2,X1] = meshgrid( -floor(res(2)/2):floor((res(2)-1)/2),...
    -floor(res(1)/2):floor((res(1)-1)/2) );