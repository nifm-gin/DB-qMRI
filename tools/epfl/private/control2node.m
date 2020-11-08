%% control2node.m
%
% Function that returns the coordinates of the node points of the curve
% (Bezier parametrization), given the control points (spline parametrization).
%
% Matthieu Guerquin-Kern, Biomedical Imaging Group / EPF Lausanne,
% 23-07-2009 (dd-mm-yyyy)

function node = control2node(control,shift)

if nargin<2
    shift = 0;
end

if shift==0 % fast computation
    node = 0.5*( circshift(control,[1 0]) + circshift(control,[2 0]) );
else
    node = control;
    H = fft(beta2((0:size(control,1)-1)+shift)).';
    node(:,1) = real(ifft(fft(control(:,1)).*H));
    node(:,2) = real(ifft(fft(control(:,2)).*H));
end