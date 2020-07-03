function [x, t] = dummyMRsignal(p, w)

if length(p) ~= length(w)
    error('Parameter vector and weight vector sizes are not the same')
end


% Constants and function
t       = (10:10:1000) *1e-3;
f       = 1e1;

func    = @(a,prop) sin(prop*2*pi*f*t) .* exp(-t/a);
%

for c = 1:length(p)
    sumfun(c,:) = func(p(c), w(c));
end
x       = abs(sum(sumfun,1));
