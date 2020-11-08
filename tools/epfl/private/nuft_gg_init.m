function st = nuft_gg_init(w, M, Msp, Mr, tau)
% This routine Initializes a structure for 2-dimensional non-uniform 
% Fourier transform using a Gaussian interpolation kernel. It implements
% Greengard's fast Gaussian gridding.
% 
% INPUT
%	w           [N,2]	angular frequencies of interest
%	M           [2]		image dimensions (M1,M2)
%	Msp         [1]		number of neighbors used in each direction
%	Mr          [2]		oversampled grid dimensions
%   tau         [2]     variance of Gaussian interpolator (optional)
% OUPUT
%	st.E1		[N]
%   st.E2x      [N]
%   st.E2y      [N]
%   st.E3x      [Msp]
%   st.E3y      [Msp]
%   st.E4       [N]   Deconvolution kernel
%
% Adapted from the NUFFT code by Will Grissom and Jeff Fessler,
% The University of Michigan, 2008.
%
% The code handles arbitrary dimensions.
%
% Copyright, Matthieu Guerquin-Kern, 2012

if numel(M)==1
    M = [M,M];
end
if numel(Mr)==1
    Mr = [Mr,Mr];
end
if any(M>Mr)
    error('the oversampled grid must be of larger dimensions!');
end
if nargin==5
    st.tau	= tau;
else % if tau is not given
    st.tau	= pi*Msp./Mr./(Mr-M/2);
end
if numel(st.tau)==1
    st.tau = [st.tau, st.tau];
end

st.M	= M;
%st.Msp	= Msp; % given by the length of E3x
%st.R   = Mr./M; % useless
st.Mr   = Mr;

% precomputations
% parts of gaussian interpolation kernel
st.E3x = exp(-((pi*(0:Msp)'/st.Mr(1)).^2)/st.tau(1)); % this is symmetric we compute only for positive coeff
st.E3y = exp(-((pi*(0:Msp)'/st.Mr(2)).^2)/st.tau(2));
%[p1,p2] = ndgrid(0:M(1)-1,0:M(2)-1);
[p1,p2] = ndgrid(-floor(M(1)/2):floor(M(1)/2+1/2)-1,-floor(M(2)/2):floor(M(2)/2+1/2)-1);
st.E4 = pi*exp(st.tau(1)*p1.^2+st.tau(2)*p2.^2)/sqrt(prod(st.tau)); % deconvolution kernel

% Nearest grid points and corresponding frequencies
st.m1 = floor(w(:,1)/2/pi*st.Mr(1));
st.m2 = floor(w(:,2)/2/pi*st.Mr(2));
w1 = 2*pi*st.m1/st.Mr(1); % nearest normalised angular frequencies on oversampled grid
w2 = 2*pi*st.m2/st.Mr(2); % nearest normalised angular frequencies on oversampled grid
st.m1 = int32(mod(st.m1,st.Mr(1))); % no reason to keep this data in double float
st.m2 = int32(mod(st.m2,st.Mr(2)));

% calculate irregular-point-dependent terms
st.E1 = exp(-((w(:,1)-w1).^2/st.tau(1) + (w(:,2)-w2).^2/st.tau(2))/4 - 1j*(w(:,1)*floor(M(1)/2)+w(:,2)*floor(M(2)/2)));
st.E2x = exp(pi*(w(:,1)-w1)./st.Mr(1)./st.tau(1));
st.E2y = exp(pi*(w(:,2)-w2)./st.Mr(2)./st.tau(2));