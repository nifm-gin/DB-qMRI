function [X, t] = ScalableSignals(p, w)

% Fabien Boux - 11/2020

narginchk(2, 2);

if size(p,2) ~= length(w)
    error('Parameter vector and weight vector sizes are not the same')
end

% Constants and function
t       = (10:10:1000) *1e-3;
f       = 50;
scalfunc = @(a,prop) sin(prop*f*t) .* exp(-t./a);

for c = 1:length(w)
    sumfun(:,:,c) = scalfunc(p(:,c), w(c));
end
X       = abs(sum(sumfun,3));
