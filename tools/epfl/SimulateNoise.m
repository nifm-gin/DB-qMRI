%% n = SimulateNoise(m, k, snr_hf, an)
%
% Generates noise to corrupt MRI data. The noise power is computed
% relatively to the average power in the highest frequencies sampled.
%
% Copyright, Matthieu Guerquin-Kern, 2012

function n = SimulateNoise(m, k, varargin)
%% Default inputs
if nargin<2
    error('SimulateNoise:TooFewInputs', ...
        'requires at least 2 inputs');
end
numvarargs = length(varargin);
if numvarargs > 5
    error('SimulateNoise:TooManyInputs', ...
        'requires at most 4 inputs');
end
optargs = {.5 false};
optargs(1:numvarargs) = varargin;
[snr_hf,an] = optargs{:};

%% Computations
k = k(:,1)+1j*k(:,2);
maxk = max([abs(real(k));abs(imag(k))]);
absk = abs(k);
ind = find((absk>0.9*maxk)&(absk<maxk));
msub = m(ind,:);
n = randn(size(m))+1j*randn(size(m));
n = n - repmat(mean(n,1),[size(m,1),1]);
n = snr_hf*std(msub(:))*n/std(n(:));

if an
    kbin = linspace(0,max(absk),15);
    for i = 2:numel(kbin)
        ind = find((abs(k(:))<kbin(i))&(abs(k(:))>=kbin(i-1)));
        dm(i-1) = std(reshape(m(ind,:),[numel(ind)*size(m,2),1]));
        dn(i-1) = std(reshape(n(ind,:),[numel(ind)*size(m,2),1]));
    end
    figure(65);h = bar(kbin(2:end),log10([dm(:),dn(:)]));title('log10 of power density');xlabel('frequency norm ||w||');legend('data','noise');
end