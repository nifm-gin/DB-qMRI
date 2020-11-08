%% DemoBrainEPI.m
%
% This script loads real scanner data from an Echo Plannar Imaging
% experiment (Cartesian k-space sampling) and precomputed receiving coil
% sensitivity maps. The full dataset is first processed and a reference
% image is reconstructed out of it. Then, a reduced dataset is loaded and
% a more challenging reconstruction takes place. The resulting image is
% compared to the reference. Note that the processing of the data includes
% a modulation such that our convention that defines the origin at the 
% upper left corner of the image is satisfied. Results using this dataset
% were presented in [1].
%
% [1] M. Guerquin-Kern, M. Häberlin, K. P. Pruessmann, and M. Unser,
% "A fast wavelet-based reconstruction method for magnetic resonance
% imaging," IEEE Transactions on Medical Imaging, vol. 30, no. 9,
% pp. 1649-1660, September 2011.
%
% Copyright, Matthieu Guerquin-Kern, 2012

clear all;close all;clc
addpath('scannerdata');

%% reconstruct full dataset
load BrainEPI_200x200_channel1_2_3_4_5_6_7_8_R1-0.mat;
load BrainEPIsens.mat;

Ncoils = size(m,2);
% Prepare Data for Recon
% getting support
support = (RSoS>3e-2*max(RSoS(:)));
support = imdilate(support,strel('disk',5));
% adjust k-space traj
k = [-k(2,:)' k(1,:)'];
% modulate data such that the origin is at the upper left corner.
m = m.*repmat(exp(1i*pi*(k(:,1)+k(:,2))),[1, Ncoils]);

[a,A] = Prepare4Recon(m, k, sensitivitites,support);
% here the preconditionner P is useless because the sensitivities have been
%estimated assuming the sum of square sensitivities is uniform.
x = a; % starting point

clear k m n w sens sensitivitites;

lambda = 1e-5*max(abs(a(:)));
ref = ReconCG(a,@(x) A(x) + lambda*x,x,30);

figure(1);imagesc(abs(ref));colormap gray;axis image;colorbar;title('image reconstructed out of the full dataset');

%% reconstruct full dataset
load BrainEPI_200x200_channel1_2_3_4_5_6_7_8_R4-33.mat;
%load BrainEPIsens.mat;

Ncoils = size(m,2);
% Prepare Data for Recon
disp('Preparing data for reconstruction');
% getting support
kernel = ones(7);kernel = kernel/sum(kernel(:));
support = (conv2(RSoS,kernel,'same')>3e-2*max(RSoS(:)));
% adjust k-space traj
k = [-k(2,:)' k(1,:)'];
% modulate data such that the origin is at the upper left corner.
m = m.*repmat(exp(1i*pi*(k(:,1)+k(:,2))),[1, Ncoils]);

[a,A] = Prepare4Recon(m, k, sensitivitites,support);
x = a; % starting point

clear k m n w sens sensitivitites;

lambda = 1e-5*max(abs(a(:)));
[x,t,d] = ReconCG(a,@(x) A(x) + lambda*x,x,30);

err = x-ref;
fprintf('\t-> Reconstruction performed with %d iterations in %.2f seconds\n',numel(t)-1,t(end));
fprintf('\t-> Reconstruction SER %.2f dB\n',-20*log10(norm(err(:))/norm(ref(:))));
figure(1);imagesc(abs(x));colormap gray;axis image;colorbar;title('reconstructed image (CG)');
figure(2);imagesc(abs(err));colormap(1-gray);axis image;colorbar;title('error map (CG) in inverted gray levels');
figure(3);semilogy(t,d,'*-');xlabel('time (s)');ylabel('residual');
