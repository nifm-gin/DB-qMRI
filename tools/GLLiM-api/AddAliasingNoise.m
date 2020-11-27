function [X_noisy, real_snr] = AddAliasingNoise(X, snr)
% Add aliasing Gaussain noise due to undersampling artifact to X according
% to SNR that is the signal-to-noise ratio, see:
% [Kara, Parameter map error due to normal noise and aliasing artifacts
% in MR fingerprinting, 2019] 
%
% Inputs: 
%   - X: N*T matrix signals to noise
%   - snr: number > 0 
%
% Output:
%   - X_noisy: N*T matrix signals = X + noise
%   - real_snr: is the effective noise
%
% Fabien Boux - 10/2020

narginchk(2, 2);

if ~any(imag(X(:)) ~= 0)
    X_noisy = X + abs(X) .*(randn(size(X)) ./snr);
    
    real_snr = max(X,[],2) ./ std(X - X_noisy, [],2);
    
else
    X_noisy = complex(real(X) + abs(real(X)) .* randn(size(X)) ./snr, ...
                      imag(X) + abs(imag(X)) .* randn(size(X)) ./snr);
                  
    real_snr = max(abs(X),[],2) ./ std(abs(X) - abs(X_noisy), [],2);
end
        