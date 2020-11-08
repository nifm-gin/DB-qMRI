%% TestGeneratingSens
%
% Computing the sensitivity maps for an array of receiving coils. Fitting the
% complex maps with two parametric models and comparing the fitting errors.
%
% Copyright, Matthieu Guerquin-Kern, 2012


%clear all;clc;
close all;
disp('Performing the test for sensitivity maps...');
Nc = 3;
mxsize = 100*[1 1];
FOV = .24;
% Computing the sensitivity maps
S = GenerateSensitivityMap(FOV, FOV./mxsize, Nc, .04, .18);
% display
figure;
for i=1:Nc
    subplot(Nc,2,(i-1)*2+1);imagesc((real(S(:,:,i))));axis image;colormap jet;colorbar;title('Bx (real part)');
    subplot(Nc,2,(i-1)*2+2);imagesc((imag(S(:,:,i))));axis image;colormap jet;colorbar;title('By (imaginary part)');
end
% Fitting the sensitivity on the SL support
DefineSL;
x = RasterizePhantom(SL,mxsize);
support = (x>1e-3*max(abs(x(:))));
support = imdilate(support,strel('disk',2));
sens_poly=cell(1,Nc);
sens_sinu=cell(1,Nc);
nrmse = zeros(2,Nc);
maxerror = zeros(2,Nc);
for i=1:Nc
	[sens_poly{i},nrmse(1,i),~,maxerror(1,i)] = SensFitting(S(:,:,i),'polynomial',8,support);
	[sens_sinus{i},nrmse(2,i),~,maxerror(2,i)] = SensFitting(S(:,:,i),'sinusoidal',8,support);
end
% Fitting accuracy on the support
fprintf('Polynomial fitting:\n * average NRMSE: %g\n * max error: %g\n',mean(nrmse(1,:)),max(maxerror(1,:)));
fprintf('Sinusoidal fitting:\n * average NRMSE: %g\n * max error: %g\n',mean(nrmse(2,:)),max(maxerror(2,:)));
