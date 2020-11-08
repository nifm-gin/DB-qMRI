%% beta2.m
% 
% Function that implements the causal B-spline of order 2.
%
% Matthieu Guerquin-Kern, Biomedical Imaging Group / EPF Lausanne,
% 23-07-2009 (dd-mm-yyyy)

function y = beta2(x)

y = zeros(size(x));

ind = find((x>0) & (x<=1));
y(ind) = x(ind).^2/2;

ind = find(x>1 & x<=2);
y(ind) = -3/2+3*x(ind)-x(ind).^2;

ind = find(x>2 & x<=3);
y(ind) = 9/2-3*x(ind)+x(ind).^2/2;