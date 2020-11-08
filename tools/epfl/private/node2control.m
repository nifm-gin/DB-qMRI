%% node2control.m
%
% Function that returns the coordinates of the control points of the curve
% (spline parametrization), given the node points (Bezier).
%
% Matthieu Guerquin-Kern, Biomedical Imaging Group / EPF Lausanne,
% 23-07-2009 (dd-mm-yyyy)

function control = node2control(node,shift)

if nargin<2
    shift = 0.5;
end

control = node;
H = fft(beta2((0:size(node,1)-1)+shift)).'; % a non null shift stabilizes the inversion
control(:,1) = real(ifft(fft(node(:,1))./H));
control(:,2) = real(ifft(fft(node(:,2))./H));