%% DemoSimuAndRecon.m
%
% This script defines a parallel MRI experiment setting, with analytically
% defined phantoms. Simulation can be performed in two different ways and
% the resulting data is corrupted by noise. This synthetic scanner data is
% prepared for reconstruction. Finally, several reconstruction methods are
% performed and the resulting reconstructions are compared to the reference
% image. Note that even in a very favorable pMRI setting, reconstruction 
% will not be perfect because of the Gibb's phenomenon which reflects the
% mismatch between the continuous nature of the MRI physics and the
% discrete nature of the model used for reconstruction [1].
%
% [1] M. Guerquin-Kern, L. Lejeune, K. P. Pruessmann, and M. Unser,
% "Realistic analytical phantoms for parallel magnetic resonance imaging,"
% IEEE Transactions on Medical Imaging, vol. 31, no. 3, pp. 626-636,
% March 2012.
%
% Copyright, Matthieu Guerquin-Kern, 2012

%simu = 'analytical';
simu = 'rasterized';


%% Defining MR setup

mxsize = 256*[1,1];
R = 4;
f_sampling = 1/1.3;
%traj = 'cartesian';
traj = 'spiral';
snr_hf = .7; % snr for the data that are outside the disc with radius 90% the highest frequency sampled
Ncoils = 1;
FOV = 0.28; % FOV width
pixelsize = FOV./mxsize;


%% Define phantom
disp('Defining analytical phantom');
%DefineSimpleBrain;
DefineBrain; % Analytical computations are pretty long with this one...
%DefineSL;Brain = SL;clear SL;
Brain.FOV = FOV*[1, 1];
ref = RasterizePhantom(Brain,mxsize);
support = (ref>5e-3*max(ref(:)));
support = imdilate(support,strel('disk',5));


%% Coils simulation
disp('Simulating coils');
sensitivity = GenerateSensitivityMap( FOV, pixelsize, Ncoils, .07, .17);


%% K-SPACE
disp('Simulating K-space traj');
switch traj
    case 'cartesian'
        w = GenerateCartesianTraj(FOV, pixelsize, f_sampling, R);
    case 'spiral'
        w = GenerateSpiralTraj( FOV, pixelsize, f_sampling, R, 10, 1, 0.031, 200, false, false);
end


%% MR SIMULATIONS
m = zeros(size(w,1),Ncoils);
switch simu
    case 'analytical'
        disp('Simulating analytical MR data');
        sens=cell(1,Ncoils);
        for indCoil = 1:Ncoils
            sens{indCoil} = SensFitting(sensitivity(:,:,indCoil),'sinusoidal',6,support);
            m(:,indCoil) = MRDataAnalytical(Brain, sens{indCoil}, w);
            [~,sensitivity(:,:,indCoil)] = RasterizePhantom(Brain,mxsize,sens{indCoil});
        end
        m = prod(mxsize)*m;
    case 'rasterized'
        disp('Simulating rasterized MR data');
        refinement = 3;
        im_rast = RasterizePhantom(Brain,refinement*mxsize);
        sens = GenerateSensitivityMap( FOV, FOV./mxsize/refinement, Ncoils, .07, .17);
        for indCoil = 1:Ncoils
            m(:,indCoil) = MRDataRasterized(sens(:,:,indCoil).*im_rast, w, FOV)/refinement^ndims(im_rast);
        end
    otherwise
        error('unknown simulation method');
end


%% Add noise
disp('Simulating noise');
n = SimulateNoise(m,w,snr_hf,true);
m = m + n;

%% Prepare Data for Recon
disp('Preparing data for reconstruction');
k = TrajInGridUnits(w, FOV, mxsize);
%figure(3);plot(k(:,1),k(:,2),'.');axis square;title('k-space traj in grid units.')
[a,A,P] = Prepare4Recon(m, k, sensitivity, support);
x = a./P./P; % starting point

%clear k m n w sens sensitivity; % cleaning data not used for reconstruction


%% Basic reconstruction
disp('Basic reconstruction');
lambda = 0*1e-3*max(abs(a(:)));
[x,t,d] = ReconCG(a, @(x) A(x) + lambda*x, a, 30, P);

err = x-ref;
fprintf('\t-> Reconstruction performed with %d iterations in %.2f seconds\n',numel(t)-1,t(end));
fprintf('\t-> Reconstruction SER %.2f dB\n',-20*log10(norm(err(:))/norm(ref(:))));
figure(1);imagesc(abs(x));colormap gray;axis image;colorbar;title('reconstructed image (CG)');
figure(2);imagesc(abs(err));colormap(1-gray);axis image;colorbar;title('error map (CG) in inverted gray levels');
figure(3);semilogy(t,d,'*-');xlabel('time (s)');ylabel('residual');
% 
% %% Non linear TV Regularized Reconstruction
% disp('Non linear TV regularized reconstruction');
% 
% [x,t,d] = ReconTV(a, A, 1.0e-3*max(abs(a(:))), a, 30, 6, P);
% 
% err = x-ref;
% fprintf('\t-> Reconstruction performed with %d iterations in %.2f seconds\n',numel(t)-1,t(end));
% fprintf('\t-> Reconstruction SER %.2f dB\n',-20*log10(norm(err(:))/norm(ref(:))));
% figure(4);imagesc(abs(x));colormap gray;axis image;colorbar;title('reconstructed image (TV)');
% figure(5);imagesc(abs(err));colormap(1-gray);axis image;colorbar;title('error map (TV) in inverted gray levels');
% figure(3);semilogy(t,d,'*-');xlabel('time (s)');ylabel('residual');
% 
% %% Non linear Wavelet Regularized Reconstruction
% disp('Non linear wavelet regularized reconstruction');
% addpath('waveletstuff/')
% 
% W = DWT(3*[1,1],mxsize, true, 'haar');
% %W = BlockDCT(8); % uncomment this if you prefer the 8x8 BlockDCT as a regularization transform
% 
% alpha = PowerIteration(A, a);
% %alpha = PowerIterationWav( @(w) W*(A(W'*(w))), W*a ); % use this one to perform FWISTA
% 
% % the alpha can depends on the MRI setting, not on the data, it can be
% % precomputed and saved for future reconstructions involving same k-space
% % trajectory, coil sensitivities and reconstruction matrix.
% 
% [x,t,d] = ReconWavFISTA(a, A, 1.5e-3*max(abs(a(:))), W, alpha, a, 100, true);
% 
% err = x-ref;
% fprintf('\t-> Reconstruction performed with %d iterations in %.2f seconds\n',numel(t)-1,t(end));
% fprintf('\t-> Reconstruction SER %.2f dB\n',-20*log10(norm(err(:))/norm(ref(:))));
% figure(6);imagesc(abs(x));colormap gray;axis image;colorbar;title('reconstructed image (wavelet)');
% figure(7);imagesc(abs(err));colormap(1-gray);axis image;colorbar;title('error map (wavelet) in inverted gray levels');
% figure(3);semilogy(t,d,'*-');xlabel('time (s)');ylabel('residual');
