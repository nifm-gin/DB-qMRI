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


%% Set some  parameters
R  = 5;    % Rank of approximation
nx = 128;  % Image size
ny = nx;   % Image size

% This addpath assumes that you are in the example folder. Otherwise make
% sure the whole toolbox is in your pathdef and forget about this line.
addpath(genpath([pwd, '/..']));

%% Create Dictionary
load('example/flip_angle_pattern.mat')
nt = length(alpha);  % Number of time frames
TR0 = 4e-3;          % Basis TR used for the TR pattern
pSSFP = 1;           % Boolean (TR = TR0, if pSSFP == 0)

% Increments of 5% for both parameters
T1 =  .3; while T1(end)<6, T1 = [T1, T1(end)*1.05]; end
T2 = .05; while T2(end)<3, T2 = [T2, T2(end)*1.05]; end

% calculate the dictionary:
D = MRF_dictionary(T1, T2, [], alpha, TR0, pSSFP, [], R);

% and add some optional parameters used for plotting in admm_recon.m
D.plot_details{1} = 'title(''T1 (s)''); caxis([0  2]); colormap parula';
D.plot_details{2} = 'title(''T2 (s)''); caxis([0 .75]); colormap parula';
D.plot_details{3} = 'title(''PD (a.u.)''); caxis([0 1.2]); colormap gray';

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

Dtmp = MRF_dictionary(.37, .13, [], alpha, TR0, pSSFP);
time_series(round(T1org*100)==10,:) = repmat(Dtmp.normalization*Dtmp.magnetization', [sum(round(T1org(:)*100)==10), 1]);

Dtmp = MRF_dictionary(1.08, .07, [], alpha, TR0, pSSFP);
time_series(round(T1org*100)==20,:) = repmat(Dtmp.normalization*Dtmp.magnetization', [sum(round(T1org(:)*100)==20), 1]);

Dtmp = MRF_dictionary(1.82, .1, [], alpha, TR0, pSSFP);
time_series(round(T1org*100)==30,:) = repmat(Dtmp.normalization*Dtmp.magnetization', [sum(round(T1org(:)*100)==30), 1]);

Dtmp = MRF_dictionary(1.4, .1, [], alpha, TR0, pSSFP);
time_series(round(T1org*100)==40,:) = repmat(Dtmp.normalization*Dtmp.magnetization', [sum(round(T1org(:)*100)==40), 1]);

Dtmp = MRF_dictionary(4.5, 2.2, [], alpha, TR0, pSSFP);
time_series(round(T1org*100)==100,:)= repmat(Dtmp.normalization*Dtmp.magnetization', [sum(round(T1org(:)*100)==100),1]);
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

figure(1); imagesc(PDorg); eval(D.plot_details{3}); title('PDorg (a.u.)'); colorbar;
figure(2); imagesc(T1org); eval(D.plot_details{1}); title('T1org (s)');    colorbar;
figure(3); imagesc(T2org); eval(D.plot_details{2}); title('T2org (s)');    colorbar;



%% Build radial trajectory with a golden angle increment
GoldenAngle = pi/((sqrt(5.)+1)/2);

% The factor sqrt(2) is contrary to the common usage in MR and is used to
% fill the corners of the Cartesian k-space in order to avoid artifacts.
kr = pi * sqrt(2) * (-1+1/nx:1/nx:1).';
phi = (0:nt-1)*GoldenAngle;

k = [];
k(:,2,:) = kr * sin(phi);
k(:,1,:) = kr * cos(phi);

%% Simulate data
E = LR_nuFFT_operator(k, [nx ny], [], [], 2); % this nuFFT operator does not incorporate a low rank transformation

% Simulate noise free data:
data = E*time_series;

% uncomment to test the reconstructions with noise of SNRin = 100 [1]:
% data = E*(time_series + (randn(size(time_series)) + 1i * randn(size(time_series)))/100);


%% Construct low rank nuFFT Operator
% The compression matrix D.u has been calculated by MRF_dictionary.m based
% on an SVD of the the dictionary. It approximates the temporal evolution
% with rank R. The default 5 nearest neighbors are used for interploation 
% with a Kaiser Bessel Kernel and oversampling with a factor of 2 is 
% employed. 
ELR = LR_nuFFT_operator(k, [nx ny], D.u, [], 2);

%% SVD Back Projection Reconstruction
% Reconstruction as proposed by 
% D. McGivney, D. Ma, H. Saybasili, Y. Jiang, and M. Griswold, ?Singular 
% Value Decomposition for Magnetic Resonance Fingerprinting in the Time 
% Domain,? IEEE Trans. Med. Imaging, vol. 33, no. 12, pp. 2311?2322, 2014.

dcomp = col(l2_norm(k,2)); % density compensation, since ELR was constructed without
dcomp(squeeze(max(abs(k(:,1,:)), abs(k(:,2,:)))) > pi) = 0; % Avoid artifacts due to the sqrt(2) factor
x = ELR' * (data .* dcomp); % Filtered back projection
x = reshape(x, [size(x,1)*size(x,2) size(x,3)]);

% Dictionary matching:
clear c idx
for q=size(x,1):-1:1
    [c(q,1),idx(q,1)] = max(x(q,:) * conj(D.magnetization), [], 2);
end
PD    = c ./ D.normalization(idx).';
PD    = reshape(PD, [nx ny]);
qMaps = D.lookup_table(idx,:);
qMaps = reshape(qMaps, [nx, ny, size(D.lookup_table,2)]);
 
figure(236); imagesc(abs(PD)); colormap gray; title('PD (a.u.)'); colorbar;
figure(237); imagesc(qMaps(:,:,1)); eval(D.plot_details{1}); title('T1_{BP} (s)');    colorbar;
figure(238); imagesc(qMaps(:,:,2)); eval(D.plot_details{2}); title('T2_{BP} (s)');    colorbar;
disp('This is the result of the SVD back projection reconstruction.');
pause;

%% Set ADMM parameters reconstruct as proposed in [1]
mu1       = 1e-3;  % ADMM penalty parameter for dictionary comparison
mu2       = 0;     % No spatial regularization for now
lambda    = 0;     % No spatial regularization for now
n_iter    = 10;    % Number of ADMM iterations
n_cg_iter = 20;    % Number of CG iterations in each ADMM iteration

% Unlike the other plots, the CG plot re-scales each individual frame
[qMaps, PD, x, r] = admm_recon(ELR, data, D, n_iter, n_cg_iter, mu1, mu2, lambda, [], 1);
return

%% Here is an example on how to add some spatial regularization
mu1       = 1e-3;   % ADMM Coupling Parameter for dictionary comparison
lambda    = 5e-5;   % Spatial regularization parameter (l21-norm)
mu2       = lambda; % ADMM Coupling Parameter for spatial regulization
n_iter    = 20;    % Number of ADMM iterations
n_cg_iter = 20;    % Number of CG iterations in each ADMM iteration

P = wavelet_operator([nx ny], 3, 'db2');
% P = finite_difference_operator([1 2]);
[qMaps, PD, x, r] = admm_recon(ELR, data, D, n_iter, n_cg_iter, mu1, mu2, lambda, P, 1);