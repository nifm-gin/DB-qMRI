%% GenerateSpiralTraj.m
%
% Generates a spiral k-space trajectory with a method adapted from [1]
% (the paper is a bit buggy).
%
% INPUTS:	* field of view in meters
%           * resolution (Nyquist distance) in meters
%           * undersampling factor along frequency encoding direction (normalized units, positive, may be lower than one)
%           * undersampling factor (normalized units, positive, may be lower than one)
%           * number of interleaves
%           * variable density factor (alpha in [1])
%           * maximum amplitude (in T/m)
%           * maximum slew rate (in T/m/s)
%           * resampling trajectory (true or false)
%           * analysis option (true or false)
%
% OUTPUT:  	* Mx2 vector of trajectory points (in rad/meters)
%
% [1]   "Simple Analytic Variable Density Spiral Design"
%       Dong-hyun Kim, Elfar Adalsteinsson, and Daniel M. Spielman
%       Magnetic Resonance in Medicine 50:214-219 (2003)
%
% Copyright, Matthieu Guerquin-Kern, 2012

function w = GenerateSpiralTraj( FOV, varargin)

%% Default parameters
if nargin==0
    FOV=.24;
end
numvarargs = length(varargin);
if numvarargs > 9
    error('GenerateSpiralTraj:TooManyInputs', ...
        'requires at most 10 inputs');
end
optargs = {FOV/128 1 1 10 4 0.031, 200 true false};
optargs(1:numvarargs) = varargin;
[res f_sampling R ninterleaves alpha gm sm rs an] = optargs{:};
if any(FOV<=res)
    error('GenerateSpiralTraj:Resolution','the field of view must be wider than the resolution');
end
if numel(res)>1
    if res(1)~=res(2)
        warning('GenerateSpiralTraj:Resolution','taking max resolution: spiral trajectories make sense for isotropic resolution');
    end
    res = max(res(:));
end
if numel(FOV)>1
    if res(1)~=res(2)
        warning('GenerateSpiralTraj:FOV','taking max FOV: spiral trajectories make sense for isotropic FOV');
    end
    FOV = max(FOV(:));
end

%% Generating first interleave
gamma = 2.678e8; % in rad/T/s
lambda = .5/res; % in m^(-1)
n = 1/(1-(1-ninterleaves*R/FOV/lambda)^(1/alpha));
w = 2*pi*n;
Tea = lambda*w/gamma/gm/(alpha+1); % in s
Tes = sqrt(lambda*w^2/sm/gamma)/(alpha/2+1); % in s
Ts2a = (Tes^((alpha+1)/(alpha/2+1))*(alpha/2+1)/Tea/(alpha+1))^(1+2/alpha); % in s

if Ts2a<Tes
    tautrans = (Ts2a/Tes).^(1/(alpha/2+1));
    tau = @(t) (t/Tes).^(1/(alpha/2+1)).*(0<=t).*(t<=Ts2a)+((t-Ts2a)/Tea + tautrans^(alpha+1)).^(1/(alpha+1)).*(t>Ts2a).*(t<=Tea).*(Tes>=Ts2a);
    Tend = Tea;
else
    tau = @(t) (t/Tes).^(1/(alpha/2+1)).*(0<=t).*(t<=Tes);
    Tend = Tes;
end
% if (w*tau(Ts2a)/alpha)^2<20
%     warning('GenerateSpiralTraj:Assumption','the assumption made for the amplitude limited case is not correctly satisfied.');
% end
k = @(t) lambda*tau(t).^alpha.*exp(1i*w*tau(t));
dt = Tea*1e-4; % in s
Dt = dt*f_sampling/FOV/abs(k(Tea)-k(Tea-dt)); % in s
t = 0:Dt:Tend; % in s
kt = k(t); % in rad

%% analysis of the gradient amplitude and slew rate
if an
    g = [0, (kt(2:end)-kt(1:end-1))/Dt/gamma];
    s = [0, (g(2:end)-g(1:end-1))/Dt];
%     figure(1);plot(tau(t),abs(g),'*-');xlabel('tau');ylabel('|G|');title('Spiral traj: gradient amplitude');
%     figure(3);plot(tau(t),abs(s),'*-');xlabel('tau');ylabel('G_x');title('Spiral traj: slew rate');
%     figure(2);plot(t,real(g),'*-');xlabel('t');ylabel('G_x');title('Spiral traj: gradient along x');
end

%% resampling (flattens sampling density for alpha=1)
if rs
    dist = [0,cumsum(abs(kt(2:end)-kt(1:end-1)))];
    kt = interp1(dist,kt,0:f_sampling/FOV:dist(end),'spline');
end

%% Generating cloned interleaves
k = kt;
for i=1:ninterleaves-1
    k = [k, kt(2:end)*exp(2*pi*1i*i/ninterleaves)];
end
w = 2*pi*[real(k(:)), imag(k(:))];

%% density analysis
if an
    kbin = linspace(0,max(abs(kt)),50);
    for i = 2:numel(kbin)
        d(i-1) = numel(find((abs(k(:))<kbin(i))&(abs(k(:))>=kbin(i-1))))/(pi*(kbin(i)^2-kbin(i-1)^2));
    end
%     figure(4);bar(kbin(2:end),d);title('sampling density');xlabel('|k|');
end