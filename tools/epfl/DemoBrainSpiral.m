%% DemoBrainSpiral
%
% This script loads real scanner data from a spiral MRI experiment. By 
% default, precomputed receiving coil sensitivity maps are loaded but the
% code to compute them can be uncommented. The full dataset is first
% processed and a reference image is reconstructed out of it. Then, a 
% reduced dataset is loaded and a more challenging reconstruction takes
% place. The resulting image is compared to the reference. Note that the
% processing of the data includes a modulation such that our convention
% that defines the origin at the upper left corner of the image is
% satisfied. Results using this dataset were presented in [1].
%
% [1] M. Guerquin-Kern, M. Häberlin, K. P. Pruessmann, and M. Unser,
% "A fast wavelet-based reconstruction method for magnetic resonance
% imaging," IEEE Transactions on Medical Imaging, vol. 30, no. 9,
% pp. 1649-1660, September 2011.
%
% Copyright, Matthieu Guerquin-Kern, 2012

clear all;close all;clc
addpath('scannerdata');

%% Get full dataset
load BrainSpiral_176_R1.mat

res = 176*[1,1];
Ncoils = size(m,2);
% Prepare Data for Recon
% adjust k-space traj
k(:,1) = -k(:,1);
% Defining center of the image
Dx = 0.43;
Dy = 0.48;
% modulate data such that the origin is at the upper left corner.
m = m.*repmat(exp(2*1i*pi*(Dx*k(:,1)+Dy*k(:,2))), [1, Ncoils]);

%% Reconstructing images for each receiving channel
% sens = ones([res Ncoils]);
% sos = zeros(res);
% for i = 1:Ncoils
%     [a,A] = Prepare4Recon(m(:,i), k, sens(:,:,i));
%     x = a; % starting point
%     lambda = 1e-5*max(abs(a(:)));
%     refcoil{i} = ReconCG(a,@(x) A(x) + lambda*x,x,100);
%     sos = sos + abs(refcoil{i}).^2;
% end
%% Estimating sensitivities
% support = sos>1e-2*max(sos(:));
% support = imdilate(support,strel('disk',4));
% sens = zeros([res, Ncoils]);
% for i = 1:Ncoils
%     sens(:,:,i) = refcoil{i}./sqrt(sos);
%     sens(:,:,i) = medfilt2(abs(sens(:,:,i)),6*[1,1]).*exp(1j*angle(sens(:,:,i))); % smooth out the sensitivities
%     figure(1);subplot(2,4,i);imagesc(abs(sens(:,:,i)));axis image;colormap jet;
% end
% save('scannerdata/BrainSpiralSens.mat','sens','support');

%% Reconstructing full dataset
load BrainSpiralSens.mat;
%support = imdilate(support,strel('disk',8)); % does not work when using
%the support. I do not know why.
[a,A] = Prepare4Recon(m(:,[1:1:8]), k, sens(:,:,[1:1:8]));
x = a;
ref = ReconCG(a,@(x) A(x)+1e-4*max(abs(a(:)))*x,x,100);
imagesc(abs(ref));axis image;colormap gray;title('reference image');
%% reconstruct undersampled dataset

load BrainSpiral_176_R4.mat
load BrainSpiralSens.mat
coilselect = 1:1:Ncoils;
% adjust k-space traj
k(:,1) = -k(:,1);
% Defining center of the image
Dx = 0.43;
Dy = 0.48;
% modulate data such that the origin is at the upper left corner.
m = m.*repmat(exp(2*1i*pi*(Dx*k(:,1)+Dy*k(:,2))), [1, Ncoils]);
[a,A] = Prepare4Recon(m(:,coilselect), k, sens(:,:,coilselect));

[x,t,d] = ReconCG(a,@(x) A(x)+1e-3*max(abs(a(:)))*x,a,100);

err = x-ref;
fprintf('\t-> Reconstruction performed with %d iterations in %.2f seconds\n',numel(t)-1,t(end));
fprintf('\t-> Reconstruction SER %.2f dB\n',-20*log10(norm(err(:))/norm(ref(:))));
figure(1);imagesc(abs(x));colormap gray;axis image;colorbar;title('reconstructed image (CG)');
figure(2);imagesc(abs(err));colormap(1-gray);axis image;colorbar;title('error map (CG) in inverted gray levels');
figure(3);semilogy(t,d,'*-');xlabel('time (s)');ylabel('residual');