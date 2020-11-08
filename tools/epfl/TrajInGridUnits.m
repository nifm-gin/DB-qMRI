%% TrajInGridUnits.m
%
% Converts a trajectory in grid units.
%
% INPUTS:	* trajectory in rad/m (Nx2 matrix)
%           * field of view in meters (can be 2x1 vector)
%           * matrix size for reconstruction, optional (can be 2x1 vector)
%
% OUTPUT:  	* trajectory in grid units (Nx2 matrix)
%           * given or suggested matrix size for reconstruction
%
% Copyright, Matthieu Guerquin-Kern, 2012

function [k,mx] = TrajInGridUnits(k, FOV, mx)

%% Default parameters

if nargin>3
    error('TrajInGridUnits:TooManyInputs','requires at most 5 inputs.');
end

if numel(FOV)==1
    FOV = FOV*[1,1];
end

%% Computations
k=k*diag(FOV)/2/pi;

if nargin<3 % trying to estimate matrix size
    modk = sqrt(k(:,1).^2+k(:,2).^2);
    xM = max(k(:,1));
    xm = min(k(:,1));
    yM = max(k(:,2));
    ym = min(k(:,2));
    x = max(abs(xM),abs(xm));
    y = max(abs(yM),abs(ym));
    if (max(modk)>1.1*max(x,y)) %  high frequency corners sampled
        disp(' we assume the trajectory is Cartesian');
        mx = round([xM-xm , yM-ym])+1;
    elseif (abs(x-y) < max(x,y)/10) % high frequency corners not sampled and square matrix size
        disp(' we assume the trajectory is radial or spiral');
        m = max(modk(:));
        if round(m)-m>=-m*eps
            mx = 2*floor(m);
        else
            mx = 2*floor(m)+1;
        end
    else
        error('TrajInGridUnits:TrajType', ...
            'could not guess the trajectory type. You will have to determine yourself desired matrix size.')
    end
end
if numel(mx)==1
    mx = mx*[1,1];
end
