%% GenerateRadialTraj.m
%
% Generates a radial k-space trajectory.
%
% INPUTS:	* field of view in meters
%           * resolution (Nyquist distance) in meters
%           * undersampling factor along frequency encoding direction (normalized units, positive, may be lower than one)
%           * undersampling factor in high frequencies(normalized units, positive, may be lower than one)
%
% OUTPUT:  	* Mx2 vector of trajectory points (in rad/meters)
%
% Copyright, Matthieu Guerquin-Kern, 2012

function w = GenerateRadialTraj( FOV, varargin)

%% Default parameters
if nargin==0
    FOV=.24;
end
numvarargs = length(varargin);
if numvarargs > 3
    error('GenerateRadialTraj:TooManyInputs', ...
        'requires at most 4 inputs');
end
optargs = {FOV/128 1 1};
optargs(1:numvarargs) = varargin;
[res f_sampling R] = optargs{:};
if any(FOV<=res)
    error('GenerateRadialTraj:Resolution','the field of view must be wider than the resolution');
end
if numel(res)>1
    warning('GenerateRadialTraj:Resolution','taking max resolution: radial trajectories make sense for isotropic resolution');
    res = max(res(:));
end
if numel(FOV)>1
    warning('GenerateRadialTraj:FOV','taking max FOV: radial trajectories make sense for isotropic FOV');
    FOV = max(FOV(:));
end

%% Generating first interleave
kt = 0:f_sampling/FOV:1/res;
kt = kt-kt(ceil((end+1)/2));

%% Generating cloned interleaves
ninterleaves = round(pi/2/atan(R/2/FOV/max(abs(kt(:)))));
k = kt;
for i=1:ninterleaves-1
    k = [k, kt*exp(pi*1i*i/ninterleaves)];
end
w = 2*pi*[real(k(:)), imag(k(:))];