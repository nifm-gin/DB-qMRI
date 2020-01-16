%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
% This a test script to demonstrate the interface of the NYU MR 
% Fingerprinting reconstruction toolbox. The test example employs a pSSFP
% [3] acquisition pattern and a radial k-space trajectory with a golden
% angle increment. 
% The two key elements of the toolbox are the low rank nuFFT operator class 
% (LR_nuFFT_operator.m) and the ADMM algorithm (admm_recon.m). The low rank
% nuFFT uses Fessler's toolbox to calculate the convolution kernel and
% multiplies it with the SVD compression matix, exploiting the commutative
% properties of the spatial FFT and the temporal low rank transformation.
% The operator acts similar to a matrix with the mtimes function (E*x) at
% its core. The function admm_recon alternately solves the constraint low
% rank inverse problem and fits SVD-series of images to the dictionary.
%
% If you use this toolbox for preparing a publication, we would highly
% apprechiate if you cite the following paper:
%
% [1] J. Asslaender, M.A. Cloos, F. Knoll, D.K. Sodickson, J.Hennig and
%     R. Lattanzi, Low Rank Alternating Direction Method of Multipliers
%     Reconstruction for MR Fingerprinting  Magn. Reson. Med., epub
%     ahead of print, 2016.
%
% [2] J. A. Fessler and B. P. Sutton, Nonuniform fast fourier
%     transforms using min-max interpolation, IEEE Trans. Signal
%     Process., vol. 51, no. 2, pp. 560-574, Feb. 2003.
%
% More details on the pSSFP sequence can be found in:
%
% [3] J. Asslaender, S. J. Glaser, and J. Hennig, ?Pseudo Steady-State Free 
%     Precession for MR-Fingerprinting,? Magn. Reson. Med., p. epub ahead 
%     of print, 2016.
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) Jakob Asslaender, August 2016
% New York University School of Medicine, Center for Biomedical Imaging
% University Medical Center Freiburg, Medical Physics
% jakob.asslaender@nyumc.org
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Simulation settings
FA      = (pi/180)* [90 ...
           10 + 50 *sin((2*pi/500) *(1:250)) + 5*randn(1,250) ...
           zeros(1,49) ...
           (10 + 50 *sin((2*pi/500) *(1:200)) + 5*randn(1,200)) /2]; % Flip angles
% TR      = 1e-3 * (10.5  + 3.5 * rand(1,500));    % Uniform between 10.5 and 14 ms
TR      = perlin(500); TR = TR(1,:);
TR      = (TR - min(TR)) ./ max(TR - min(TR));
TR      = 1e-3 * (10.5  + 3.5 * TR);

FA = repmat(FA,1,2);
TR = repmat(TR,1,2);

FA = FA(1:850);
TR = TR(1:850);


%% Set some  parameters
R  = [];    % Rank of approximation
nx = 128;  % Image size
ny = nx;   % Image size

% Adding to path
addpath(genpath('functions'))
addpath(genpath('tools'))

%% Create Dictionaries
load('tools/NYU_MRF_Recon/example/flip_angle_pattern.mat')
nt = length(alpha);  % Number of time frames
TR0 = 4e-3;          % Basis TR used for the TR pattern
pSSFP = 1;           % Boolean (TR = TR0, if pSSFP == 0)

% Increments of 5% for both parameters
T1  =  .3; while T1(end)<6, T1 = [T1, T1(end)*1.09]; end
T2  = .05; while length(T2)<length(T1), T2 = [T2, T2(end)*1.12]; end
T1  = repmat(T1, 1, length(T2));
T2  = repelem(T2, 1, length(T2));

% calculate the grid dictionary:
Dgrid = MRF_dictionary(T1, T2, [], FA, TR);%, pSSFP, [], R);


T  	= net(scramble(sobolset(2),'MatousekAffineOwen'),length(T1))';
T1  = min(T1) + (max(T1) - min(T1)) * T(1,:);
T2  = min(T2) + (max(T2) - min(T2)) * T(2,:);

% calculate the quasi-random dictionary:
Dqmc = MRF_dictionary(T1, T2, [], FA, TR);%, TR0, pSSFP, [], R);

% and add some optional parameters used for plotting in admm_recon.m
% Dgrid.plot_details{1} = 'title(''T1 (s)''); caxis([0  2]); colormap parula';
% Dgrid.plot_details{2} = 'title(''T2 (s)''); caxis([0 .75]); colormap parula';
% Dgrid.plot_details{3} = 'title(''PD (a.u.)''); caxis([0 1.2]); colormap gray';


%% Create numerical phantom
% This numerical phantom is based on the Shepp-Logan phantom and relaxation
% paramters around the values of brain tissue at 3T are used. Please note,
% that a Shepp Logan phantom does not provide a realistic scenario, if
% spatial regularization is applied and we use it in this demonstration for
% the sake of convenience.

T1org = phantom(nx,ny);
T2org = T1org;
PDorg = T1org;
time_series = repmat(T1org(:), [1,nt]);

Dtmp = MRF_dictionary(.37, .13, [], FA, TR);
time_series(round(T1org*100)==10,:) = repmat(Dtmp.normalization*Dtmp.magnetization', [sum(round(T1org(:)*100)==10), 1]);
% time_series(round(T1org*100)==10,:) = repmat(Dtmp.magnetization', [sum(round(T1org(:)*100)==10), 1]);

Dtmp = MRF_dictionary(1.08, .07, [], FA, TR);
time_series(round(T1org*100)==20,:) = repmat(Dtmp.normalization*Dtmp.magnetization', [sum(round(T1org(:)*100)==20), 1]);
% time_series(round(T1org*100)==20,:) = repmat(Dtmp.magnetization', [sum(round(T1org(:)*100)==20), 1]);

Dtmp = MRF_dictionary(1.82, .1, [], FA, TR);
time_series(round(T1org*100)==30,:) = repmat(Dtmp.normalization*Dtmp.magnetization', [sum(round(T1org(:)*100)==30), 1]);
% time_series(round(T1org*100)==30,:) = repmat(Dtmp.magnetization', [sum(round(T1org(:)*100)==30), 1]);

Dtmp = MRF_dictionary(1.4, .1, [], FA, TR);
time_series(round(T1org*100)==40,:) = repmat(Dtmp.normalization*Dtmp.magnetization', [sum(round(T1org(:)*100)==40), 1]);
% time_series(round(T1org*100)==40,:) = repmat(Dtmp.magnetization', [sum(round(T1org(:)*100)==40), 1]);

Dtmp = MRF_dictionary(4.5, 2.2, [], FA, TR);
time_series(round(T1org*100)==100,:)= repmat(Dtmp.normalization*Dtmp.magnetization', [sum(round(T1org(:)*100)==100),1]);
% time_series(round(T1org*100)==100,:)= repmat(Dtmp.magnetization', [sum(round(T1org(:)*100)==100),1]);
time_series = reshape(time_series, nx,ny,nt);

PDorg(PDorg>0) = 1;

T1org(round(T1org*100)==10) = .37; %s
T1org(round(T1org*100)==20) = 1.08; %s
T1org(round(T1org*100)==30) = 1.82; %s
T1org(round(T1org*100)==40) = 1.4; %s
T1org(round(T1org*100)==100) =  4.5; %s

T2org(round(T2org*100)==10) = .13; %s
T2org(round(T2org*100)==20) = .07; %s
T2org(round(T2org*100)==30) = .1; %s
T2org(round(T2org*100)==40) = .1; %s
T2org(round(T2org*100) ==100) =  2.2; %s

% figure(1); imagesc(PDorg); eval(D.plot_details{3}); title('PDorg (a.u.)'); colorbar;
% figure(2); imagesc(T1org); eval(D.plot_details{1}); title('T1org (s)');    colorbar;
% figure(3); imagesc(T2org); eval(D.plot_details{2}); title('T2org (s)');    colorbar;


%% Estimates

DicoG{1}.MRSignals = abs((Dgrid.normalization.*Dgrid.magnetization)'); 
% DicoG{1}.MRSignals = abs(Dgrid.magnetization'); 
DicoG{1}.Parameters.Par = Dgrid.lookup_table;

DicoR{1}.MRSignals = AddNoise(abs((Dqmc.normalization.*Dqmc.magnetization)'), 60); 
% DicoR{1}.MRSignals = AddNoise( abs(Dqmc.magnetization'), 60); 
DicoR{1}.Parameters.Par = Dqmc.lookup_table;

Parameters.K = 50;
Parameters.cstr.Sigma  = 'd*';
Parameters.cstr.Gammat = ''; 
Parameters.cstr.Gammaw = '';
Parameters.Lw = 0;

Xtest = AddNoise(abs(time_series), 60);

% Perform DBM
Estim_dbm   = AnalyzeMRImages(Xtest,DicoG,'DBM',[],[]);

% mNRMSE_grid(snr,f,:) = Estim.GridSearch.Errors.Nrmse;
% mRMSE_grid(snr,f,:)	= Estim.GridSearch.Errors.Rmse;
% mMAE_grid(snr,f,:) 	= Estim.GridSearch.Errors.Mae;

% Perform DBL
Estim_dbl   = AnalyzeMRImages(Xtest,DicoR,'DBL',Parameters,[]);

% mNRMSE_gllim(snr,f,:) = Estim.Regression.Errors.Nrmse;
% mRMSE_gllim(snr,f,:) = Estim.Regression.Errors.Rmse;
% mMAE_gllim(snr,f,:) = Estim.Regression.Errors.Mae;


%%
mx = max(Dgrid.lookup_table);
mn = min(Dgrid.lookup_table);
bds{1} = [mn(1) mx(1)];
bds{2} = [mn(2) mx(2)];

bds_err = [0 3e-2];

mask = logical(T1org == 0 + T2org == 0);

figure
subplot(2,5,1)
imagesc(T1org .* mask, bds{1})
axis off; axis image; colorbar
title('Real')
ylabel('T_1')

subplot(2,5,6)
imagesc(T2org .* mask, bds{2})
axis off; axis image; colorbar
ylabel('T_2')

Tt{1} = T1org;
Tt{2} = T2org;
for i = 1:2
    subplot(2,5,(i-1)*5+2)
    imagesc(Estim_dbm.GridSearch.Y(:,:,i) .* mask, bds{i})
    axis off; axis image; colorbar
    if i ==1, title('DBM'); end
    
    subplot(2,5,(i-1)*5+3)
    imagesc((Tt{i} - Estim_dbm.GridSearch.Y(:,:,i)).^2 .* mask, bds_err)
    axis off; axis image; colorbar
    title(['RMSE = ' num2str(nanmean(reshape((Tt{i} - Estim_dbm.GridSearch.Y(:,:,i)).^2 .* mask,1,[])).^.5)]) 
    
    subplot(2,5,(i-1)*5+4)
    imagesc(Estim_dbl.Regression.Y(:,:,i) .* mask, bds{i})
    axis off; axis image; colorbar
    if i ==1, title('DBL'); end
    
    subplot(2,5,(i-1)*5+5)
    imagesc((Tt{i} - Estim_dbl.Regression.Y(:,:,i)).^2 .* mask, bds_err)
    axis off; axis image; colorbar
    title(['RMSE = ' num2str(nanmean(reshape((Tt{i} - Estim_dbl.Regression.Y(:,:,i)).^2 .* mask,1,[])).^.5)]) 
end

