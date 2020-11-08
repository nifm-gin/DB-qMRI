%% SensFitting.m
%
% Function that fits a given sensitivity profile by a 2D polynomial of
% given degree or a linear combination of complex sinusoids.
%
% INPUT:    * sens: sensitivity map
%           * model: 'polynomial' or 'sinusoidal'
%           * param: parameter of the model (polynomial degree or bandwith)
%           * support (optional): mask for the sensitivity fitting
%
% OUTPUT:   * s: structure defining the continuous sensitivity profile
%           * nrmse: normalized root mean square error
%           * ser: signal to error ratio of the approximation
%           * maxerror: maximal error of approximation inside the support
%           * condition number of the matrix to be inverted (gives an
%               indication of the accuracy of the results)
%
% See: subfunction Fitting
%
% Matthieu Guerquin-Kern, Biomedical Imaging Group / EPF Lausanne,
% 19-02-2011 (dd-mm-yyyy)

function [s,nrmse,ser,maxerror,condi] = SensFitting(sens,model,param,support)

res = size(sens);
if nargin<4
    support = ones(res);
end

ind = find(support);    % indexes corresponding to the support

% Define the points inside the support
[X1,X2] = GenerateFullCart2DGrid(res);
R = [X1(ind)'/res(1);X2(ind)'/res(2)];
clear X1 X2;

% Matrix related to the points of the support
switch model
    case 'polynomial'
        M = Polynomial2DMatrix(R,param);
    case 'sinusoidal'
        M = Sinusoidal2DMatrix(R,param);
end

[coeff,nrmse,ser,maxerror,condi] = Fitting(M,sens(ind),'accurate');

s.model = model;
s.data = coeff;
s.param = param;

%% Fitting.m
%
% Function that fits the parameters for a given linear model.
%
% INPUT:    * M:        matrix that maps the coefficients to the data
%           * data:     data points to be fitted
%           * method:   'fast' or 'accurate' (default)
%
% OUTPUT:   * s: (N,1) vector of fitted coefficients
%           * nrmse: normalized root mean square error
%           * ser: signal to error ratio of the approximation
%           * maxerror: maximal error of approximation inside the support
%           * condition number of the matrix to be inverted (gives an
%               indication of the accuracy of the results)
%
% See: SensFitting
%
% Matthieu Guerquin-Kern, Biomedical Imaging Group / EPF Lausanne,
% 31-01-2011 (dd-mm-yyyy)

function [s,nrmse,ser,maxerror,condi] = Fitting(M,data,method)

if nargin<3
    method = 'accurate'; % default method
end

% Compute the parameters
switch method
    case 'accurate'
        s = M\data; % linear system to invert
    case 'fast'
        H = M'*M;
        s = M'*data;
        s = H\s;
    otherwise
        error('fitting method must be either ''fast'' or ''accurate''.');
end

if nargout>1
    residue = data-M*s;
    mse = residue'*residue/numel(residue);
    nrmse = sqrt(mse)/(max(abs(data)-min(abs(data))));
    ser = (data'*data)/mse/numel(residue);
    maxerror = max(abs(residue))/(max(abs(data)-min(abs(data))));
    switch method
        case 'accurate'
            condi = cond(M,2);
        case 'fast'
            condi = cond(H,2);
    end
end
