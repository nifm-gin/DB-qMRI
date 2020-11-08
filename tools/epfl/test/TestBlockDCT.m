%% TestBlockDCT.m
%
% Checking the implementation of the BlockDCT class. It mimics the block
% discrete cosine transform matrices thanks to overloaded functions.
%
% Copyright, Matthieu Guerquin-Kern, 2012

disp('Performing the test for the BlockDCT class...');
close all;
mxsize = 256*[1 1];
DefineBrain;
x = RasterizePhantom(Brain,mxsize);
T = BlockDCT(8);
c = T*x;
z = inv(T)*c;
t = T'*c;
disp('the following values should be close to 0:');
fprintf('\tchecking Parseval''s theorem: %g\n',abs(norm(c(:))-norm(x(:)))/norm(x(:))); % according to Parceval's theorem, if the transform unitary, this should be 0
fprintf('\tchecking accuracy 1: %g\n',norm(x(:)-z(:))/norm(x(:)));
fprintf('\tchecking accuracy 2: %g\n',norm(x(:)-t(:))/norm(x(:)));

figure(1);imagesc(x);axis image;colormap gray;title('original image');colorbar;
figure(2);imagesc(c);axis image;colormap jet;title('DCT coefficients');colorbar;
%figure(3);imagesc(abs(x-z));axis image;colormap(1-gray);title('error');colorbar;
