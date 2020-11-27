function [X_noisy, real_snr] = AddNoise(X, snr)
% Add gaussian white noise to X according to SNR that is the
% signal-to-noise ratio
%
% Inputs: 
%   - X: N*T matrix signals to noise
%   - snr: number > 0 
%
% Output:
%   - X_noisy: N*T matrix signals = X + noise
%   - real_snr: is the effective noise
%
% Fabien Boux - 12/2017

narginchk(2, 2);

% real signals
if ~any(imag(X(:)) ~= 0)
    X_noisy = abs(X + randn(size(X)) .* repmat(max(abs(X),[],2) ./ snr, 1,size(X,2)));
    
    real_snr = max(X,[],2) ./ std(X - X_noisy, [],2);
    
% complex signals
else
    X_noisy = complex(real(X) + randn(size(X)) .* max(abs(real(X)),[],2) ./ snr, ...
                      imag(X) + randn(size(X)) .* max(abs(imag(X)),[],2) ./ snr);
                  
    real_snr = max(abs(X),[],2) ./ std(abs(X) - abs(X_noisy), [],2);
end
        