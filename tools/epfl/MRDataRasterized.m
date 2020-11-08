%% MRDataRasterized.m
%
% Function that computes kspace data given the rasterized
% spatial data.
%
% To reduce aliasing artifacts, the resolution of the spatial data should
% be as large as possible.
%
% Copyright, Matthieu Guerquin-Kern, 2012

function m = MRDataRasterized(x,w,FOV)
if numel(FOV)==1
    FOV = FOV*[1 1];
end
param.k = w*diag(FOV)/2/pi;
param.res = size(x);
param.method = 'nuft gg';
m = E0(x,param);