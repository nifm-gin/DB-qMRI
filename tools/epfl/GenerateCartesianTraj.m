%% GenerateCartesianTraj.m
%
% Generates a Cartesian k-space trajectory.
%
%
% INPUTS:	* field of view in meters (can be 2x1 vector)
%           * resolution (Nyquist distance) in meters (can be 2x1 vector)
%           * undersampling factor along frequency encoding direction (normalized units, positive, may be lower than one)
%           * undersampling factor along phase encoding direction (normalized units, positive, may be lower than one)
%
% OUTPUT:  	* Mx2 vector of trajectory points (in rad/meters)
%
% Copyright, Matthieu Guerquin-Kern, 2012

function w = GenerateCartesianTraj(FOV, varargin)

%% Default parameters
if nargin==0
    FOV=.24;
end
numvarargs = length(varargin);
if numvarargs > 3
    error('GenerateCartesianTraj:TooManyInputs', ...
        'requires at most 4 inputs');
end
optargs = {FOV/128 1 1};
optargs(1:numvarargs) = varargin;
[res f_sampling R] = optargs{:};
if any(FOV<=res)
    error('GenerateCartesianTraj:Resolution','the field of view must be wider than the resolution');
end
if numel(res)==1
    res = res*[1,1];
end
if numel(FOV)==1
    FOV = FOV*[1,1];
end

%% Generating trajectory

x = 0:f_sampling/FOV(1):1/res(1)-f_sampling/FOV(1);
x = x-x(ceil((end+1)/2));

y = 0:R/FOV(2):1/res(2)-R/FOV(2);
y = y-y(ceil((end+1)/2));

w = zeros(numel(x)*numel(y),2);
[wx, wy] = ndgrid(2*pi*x,2*pi*y);
w(:,1) = wx(:);
w(:,2) = wy(:);
