%% [x, t, delta] = ReconWavFISTA(a, A, mu, W, alpha, x0, Nit, RS)
%
% Performs wavelet-regularized reconstruction [1] trying to solve the linear
% system
%                       Ax = a,
% with A a Hermitian symmetric matrix.
%
% The DWT object W provided for regularization should perform an
% orthogonal transform.
%
% See the classes DWT and WAVELET.
%
% [1] M. Guerquin-Kern, M. Häberlin, K. P. Pruessmann, and M. Unser,
% "A fast wavelet-based reconstruction method for magnetic resonance
% imaging," IEEE Transactions on Medical Imaging, vol. 30, no. 9,
% pp. 1649-1660, September 2011.
%
% Copyright, Matthieu Guerquin-Kern, 2012

function [x,t,delta] = ReconWavFISTA(a, A, mu, varargin)
%% Default inputs
if nargin<3
    error('ReconWavFISTA:TooFewInputs', ...
        'requires at least 3 inputs');
end
numvarargs = length(varargin);
if numvarargs > 5
    error('ReconWavFISTA:TooManyInputs', ...
        'requires at most 8 inputs');
end
optargs = {DWT(3*[1,1], size(a), true, 'haar'), 1, a, 100, true};
optargs(1:numvarargs) = varargin;
[W, alpha, x0, Nit, RS] = optargs{:};

%% Initializations
t = zeros(Nit+1,1);
delta = ones(Nit+1,1);
h = waitbar(0,'FISTA wavelet-based reconstruction...');

w = W*x0;
x = W'*w;
Ax = A(x);
Ay = Ax;
wy = w;

tau = .9 ./alpha; % cf. Combette
if length(tau)>1&&length(tau)~=length(mu)
    mu = mu(1)*ones(size(tau));mu(end) = 0;
end
WTinv = inv(W');

s = 1;

%% Computations
if RS % perform FISTA+RS
    STEP = 1;
    STEP_MAX = 20;
    xold = x;
    Axold = Ax;
    shift = round((2.^3)*(rand(Nit,2)-0.5));
    Cost = @(x,Ax,muw) real(Ax(:)'*x(:)) - 2*real(x(:)'*a(:)) + sum(abs(muw(:)));
    Cx = Cost(x,Ax,mu*w);
    for i=1:Nit
        t0 = clock();
        w = WTinv*(circshift(x,shift(i,:)));
        dMFISTA = WTinv*(circshift(x-xold,shift(i,:)));
        sold = s;
        s = (1+sqrt(1+4*sold^2))/2;
        kappa = (STEP<STEP_MAX)*min(1,(sold-1)/s);
        wy = w + kappa*dMFISTA;
        Ay = (1+kappa)*Ax - kappa*Axold;
        w = SoftThresh(wy + tau*(W*(circshift(a-Ay,shift(i,:)))),mu.*tau);
        dISTA = w - wy;
        xold = x;
        x = circshift(W'*w,-shift(i,:));
        Axold = Ax;
        Ax = A(x);
        Cxold = Cx;
        Cx = Cost(x,Ax,mu*w);
        if (Cxold<Cx)
            STEP = STEP+1;
            %disp(['increasing... iteration ' num2str(i) ' STEP:' num2str(STEP)]);
        end
        if STEP==STEP_MAX
            tau = 1*min(tau)*ones(size(tau));
            STEP = STEP+1;
        end
        t(i+1) = t(i) + etime(clock(),t0);
        delta(i+1) = norm(dISTA(:),2)/norm(w(:),2);
        waitbar(i/Nit,h);
    end
else % perform FISTA
    for i=1:Nit
        t0 = clock();
        dFISTA = SoftThresh(wy + tau*(W*(a-Ay)),mu.*tau) - w;
        
        w = w + dFISTA;
        dISTA = w - wy;
        % prepare future steps
        x = W'*(w);
        Axold = Ax;
        Ax = A(x);
        
        sold = s;
        s = (1+sqrt(1+4*sold^2))/2;
        kappa = min(1,(sold-1)/s);
        
        wy = w + kappa*dFISTA;
        Ay = (1+kappa)*Ax - kappa*Axold;
        t(i+1) = t(i) + etime(clock(),t0);
        delta(i+1) = norm(dISTA(:),2)/norm(w(:),2);
        waitbar(i/Nit,h);
    end
end
close(h);
