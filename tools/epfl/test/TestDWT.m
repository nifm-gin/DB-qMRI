%% TestDWT.m
%
% Checking the implementation of the DWT and WAVELET classes. They mimic
% discrete wavelet transform matrices and, respectively, vectors of wavelet
% coefficients thanks to overloaded functions.
%
% Copyright, Matthieu Guerquin-Kern, 2012

close all;
addpath('waveletstuff/');
disp('Performing the test for the DWT class...');
mxsize = 256*[1 1];
DefineBrain;
x = RasterizePhantom(Brain,mxsize);
W = DWT(3*[1,1],mxsize, true, 'haar'); % this wavelet transform is orthonormal
w = W*x;
z = inv(W)*w;
t = W'*w;
disp('the following values should be close to 0:');
fprintf('\tchecking Parseval''s theorem: %g\n',abs(norm(w(:))-norm(x(:)))/norm(x(:))); % according to Parceval's theorem, if the transform is orthonormal, this should be 0 
fprintf('\tchecking accuracy 1: %g\n',norm(x(:)-z(:))/norm(x(:)));
fprintf('\tchecking accuracy 2: %g\n',norm(x(:)-t(:))/norm(x(:)));

figure(1);imagesc(x);axis image;colormap gray;title('original image');colorbar;
figure(2);imagesc(w);
% figure(3);imagesc(abs(x-z));axis image;colormap(1-gray);title('error');colorbar;
