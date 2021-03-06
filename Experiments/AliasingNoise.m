
%% Description
%
% We investigate the impact of aliasing noise due to undersampling 
% artifacts
%
% Fabien Boux - 11/2020

Init
disp('Running experiment AliasingNoise.m')


%% Settings

display = 0; %for intermediate figures
backup  = 1;

nb_param = 3;   %number of parameters: 2 for T1 and T2 estimates, 
                % and 3 for T1, T2 and Df estimates
nb_signals = [16^3 61^3]; %correspond to: [4096 226981]
S       = 1000;

%model settings
Model_.K = 100;
%one can specify a file of previous saved models, else empty
load_model = [];


%% Create phantom

FOV     = 0.28; %field of view (m)
res     = 128;  %resolution

%generate mask
[ref, support, mask] = define_phantom(FOV, res);

if display == 1
    figure
    subplot(231)
    imagesc(ref)
    colormap jet; axis image; axis off
    title('Labels')
end

%define Gaussian means for parameter values
T1 = nan(size(ref)); T1_ = (1000:400:3000) *1e-3; %(s)
T2 = nan(size(ref)); T2_ = (60:60:1000) *1e-3; %(s)
T1 = nan(size(ref)); T1_ = (1500:300:3000) *1e-3; %(s)
T2 = nan(size(ref)); T2_ = (60:60:1000) *1e-3; %(s)
df = nan(size(ref));

%generate T1 and T2 values
u = unique(ref);
u(u == 0) = [];
for r = 1:numel(u)
    T1(ref == u(r)) = T1_(r);
    T2(ref == u(r)) = T2_(r);
end
T1 = T1 .* (1 + 1/8 * rand(size(T1)));
T2 = T2 .* (1 + 1/2 * rand(size(T2)));


%generate off-resonance values
max_df = 100 *pi/180; %(rad)
df_ = -max_df:(2*max_df)/size(df,1):max_df;
for i = 1:size(df,1)
    for j = 1:size(df,2)
        df(i,j) = df_(i) + df_(j);
    end
end

if display == 1
    subplot(234)
    imagesc(T1)
    colormap gray; axis image; axis off
    title('T1')
    subplot(235)
    imagesc(T2)
    colormap gray; axis image; axis off
    title('T2')
    subplot(236)
    imagesc(df)
    colormap gray; axis image; axis off
    title('df')
end


%% Simulate signals

% From Dan Ma's paper:

%simulation settings
sd      = 2; %standard deviation for FA
FA      = (pi/180)* [90 ... %Flip anglas (rad)
           10 + 50 *sin((2*pi/500) *(1:250)) + sd*randn(1,250) ...
           zeros(1,50) ...
           (10 + 50 *sin((2*pi/500) *(1:250)) + sd*randn(1,250)) /2 ...
           zeros(1,50) ...
           10 + 50 *sin((2*pi/500) *(1:250)) + sd*randn(1,250) ...
           zeros(1,50) ...
           (10 + 50 *sin((2*pi/500) *(1:99)) + sd*randn(1,99)) /2]; % Flip angles

TR      = perlin(length(FA)+100); TR = TR(1,:); % Repetition Time (sec)
TR      = filter(gausswin(50),1,TR);
TR      = TR(51:1050);
TR  =    1e-3 * (10.5 + 3.5* (TR - min(TR)) ./ max(TR - min(TR)));

TR      = TR(1:S);
FA      = FA(1:S);
%load('inputs/patterns.mat', 'FA','TR')

%img settings
N       = size(ref,1);

%simulation
if nb_param == 2
    Acq     = MRF_dictionary(reshape(T1,[],1), reshape(T2,[],1), [], FA, TR); 
elseif nb_param == 3
    Acq     = MRF_dictionary(reshape(T1,[],1), reshape(T2,[],1), reshape(df,[],1), FA, TR); 
end
Xacq    = (Acq.normalization.*Acq.magnetization).';
Iref    = reshape(Xacq, N,N, length(TR));

Xacq(any(isnan(Xacq)'),:) = 0;

% no thermal noise added
% Xacq    = AddNoise(Xacq, snr);

%generate images
Iacq    = reshape(Xacq, N,N, length(TR));

if display == 1
    figure
    subplot(211)
    plot(FA *180/pi); ylabel('FA (degrees)')
    subplot(212)
    plot(TR*1e3); ylabel('TR (ms)')
    xlabel('Repetitions')
end


%% Generate trajectories

% From Dan Ma's paper:
% The variable density spiral-out trajectory was designed to have 5.8 ms
% readout time in each TR and to have zero and first moment gradient
% compensation using minimum-time gradient design49. This trajectory
% required one interleaf to sample the inner 10×10 region, while 48
% interleaves were required to fully sample the outer portions of k-space.
% During acquisition, the spiral trajectory rotated 7.5° from one time
% point to the next, so that each time point had a slightly different
% spatial encoding

interleaves = [1 2:2:50];


%% Dico generation & Regression training

int_T1  = 1e-3 * [200 3000];	% between 20 and 3000 ms
int_T2  = 1e-3 * [20  300];    % between 20 and 3000 ms
int_df  = pi/180 * [-200 200];  % +/- 400 Hz


for d = 1:numel(nb_signals)
    
    clear X Y
    nb_step = round(nb_signals(d)^(1/nb_param));
    v1  = int_T1(1) : (int_T1(2) - int_T1(1)) / (nb_step-1) : int_T1(2);
    v2  = int_T2(1) : (int_T2(2) - int_T2(1)) / (nb_step-1) : int_T2(2);
    v3  = int_df(1) : (int_df(2) - int_df(1)) / (nb_step-1) : int_df(2);
    if nb_param == 2
        Y(:,1)  = repmat(v1, 1, length(v2));
        Y(:,2)  = repelem(v2, 1, length(v2));
        D       = MRF_dictionary(Y(:,1), Y(:,2), [], FA, TR); 
    elseif nb_param == 3
        Y(:,1)  = repmat(v1, 1, length(v2)*length(v3));
        Y(:,2)  = repmat(repelem(v2, 1, length(v1)), 1,length(v3));
        Y(:,3)  = repelem(v3, 1, length(v2)*length(v3));
        D       = MRF_dictionary(Y(:,1), Y(:,2), Y(:,3), FA, TR); 
    end
    X       = (D.normalization.*D.magnetization).';
    Dico_DBM{d} = FormatDico(X, Y(:,1:nb_param));
    disp('Dictionary for matching generated')
    
    if isempty(load_model)
        % Compute training dataset
        clear X Y
        Y       = net(scramble(sobolset(nb_param),'MatousekAffineOwen'),nb_signals(d));
        Y(:,1)  = int_T1(1) + (int_T1(2) - int_T1(1)) * Y(:,1);
        Y(:,2)  = int_T2(1) + (int_T2(2) - int_T2(1)) * Y(:,2); 
        if nb_param == 2
            D       = MRF_dictionary(Y(:,1), Y(:,2), [], FA, TR); 
        elseif nb_param == 3
            Y(:,3)  = int_df(1) + (int_df(2) - int_df(1)) * Y(:,3);
            D       = MRF_dictionary(Y(:,1), Y(:,2), Y(:,3), FA, TR); 
        end 
        X       = (D.normalization.*D.magnetization).';
        Dico_DBL = FormatDico([real(X) imag(X)], Y(:,1:nb_param));
        disp('Dictionary for learning generated')

        if size(Dico_DBM{d}.MRSignals,1) ~= size(Dico_DBL.MRSignals,1)
            warning('Sizes are not equals')
        end

        % Learn models
        [~, Model{d}] = AnalyzeMRImages([],Dico_DBL,'DB-SL',Model_);
        disp('DB-SL model learnt')

        [~, NeuralNet{d}] = AnalyzeMRImages([],Dico_DBL,'DB-DL',Model_);
        disp('DB-DL model learnt')
    end
end


%%

c = 0;
for nb_interleaves = interleaves 

    disp(['Subsampling factor = ' num2str(nb_interleaves)])
    
    SNR     = 1./(3*nb_interleaves-2)*100;
    [Irec,real_snr] = AddAliasingNoise(reshape(Iacq,[],size(Iacq,3)), SNR);
    Irec    = reshape(Irec, size(Iacq,1),size(Iacq,2),size(Iacq,3));
    
    if S >= 430
        imgrec{nb_interleaves,1} = Irec(:,:,5);
        imgrec{nb_interleaves,2} = Irec(:,:,140);
        imgrec{nb_interleaves,3} = Irec(:,:,430);
        sigrec{nb_interleaves,1} = squeeze(Irec(61,31,:));
        sigrec{nb_interleaves,2} = squeeze(Irec(52,37,:));
        sigrec{nb_interleaves,3} = squeeze(Irec(65,51,:));
    end
    
    if display == 1
        figure
        subplot(321)
        hold on
        plot(abs(squeeze(Irec(round(res/2),round(res/2*1.1),:))'))
        plot(abs(squeeze(Iacq(round(res/2),round(res/2),:))'))

        subplot(322)
        hold on
        plot(abs(squeeze(Irec(100,120,:))'))
        plot(abs(squeeze(Iacq(100,120,:))'))

        bds = [0 0.2];

        i = randi(S);
        subplot(334)
        imagesc(abs(Iref(:,:,i)), bds); colormap gray; axis image; axis off
        subplot(335)
        imagesc(abs(Iacq(:,:,i)), bds); colormap gray; axis image; axis off
        subplot(336)
        imagesc(abs(Irec(:,:,i))); colormap gray; axis image; axis off

        i = randi(S);
        subplot(337)
        imagesc(abs(Iref(:,:,i)), bds); colormap gray; axis image; axis off
        subplot(338)
        imagesc(abs(Iacq(:,:,i)), bds); colormap gray; axis image; axis off
        subplot(339)
        imagesc(abs(Irec(:,:,i))); colormap gray; axis image; axis off
    end
    clear D
    
    % Generate test data
    if nb_param == 2
        Ytest 	= [reshape(T1,[],1), reshape(T2,[],1)];
    elseif nb_param == 3
        Ytest 	= [reshape(T1,[],1), reshape(T2,[],1), reshape(df,[],1)];
    end
    Xtest   = reshape(Irec,[],size(Irec,3));
    
    nan_val = any(isnan(Ytest)');
    Ytest(nan_val,:) = [];
    Xtest(nan_val,:) = [];
    
	c = c + 1;

    for d = 1:numel(nb_signals)
        
        % Perform DBM
        Dic{1} = Dico_DBM{d};
        Estim   = AnalyzeMRImages(Xtest,Dic,'DBM',[],Ytest(:,1:nb_param));

        RMSE_grid(c,:,d)	= Estim.GridSearch.Errors.Rmse;
        MAE_grid(c,:,d) 	= Estim.GridSearch.Errors.Mae;

        % Perform DBSL
        Estim   = AnalyzeMRImages([real(Xtest) imag(Xtest)],[],'DB-SL',Model{d},Ytest(:,1:nb_param));

        RMSE_gllim(c,:,d)  = Estim.Regression.Errors.Rmse;
        MAE_gllim(c,:,d)   = Estim.Regression.Errors.Mae;
        
        % Perform DBDL
        Estim   = AnalyzeMRImages([real(Xtest) imag(Xtest)],[],'DB-DL',NeuralNet{d},Ytest(:,1:nb_param));

        RMSE_nn(c,:,d)      = Estim.Regression.Errors.Rmse;
        MAE_nn(c,:,d)       = Estim.Regression.Errors.Mae;
    end
end


%% Backup

if backup == 1
    clear X* Y* Dico* Acq Traj 
    save(['temp/' 'AlisasingNoise'])
end


%%
    
colors = [          0    0.4470    0.7410
            0.8500    0.3250    0.0980
            0.9290    0.6940    0.1250
            0.4940    0.1840    0.5560
            0.4660    0.6740    0.1880
            0.3010    0.7450    0.9330
            0.6350    0.0780    0.1840];

fig = figure;
subplot(545)
plot(FA *180/pi); ylabel('FA (degrees)')
ylim([0 90])
subplot(549)
plot(TR*1e3); ylabel('TR (ms)')
ylim([10 14.5])
xlabel('Repetitions') 


sig = 2;
subsampling_factors = [8 16 48];

if S >= 430
    img = 1;
    bds = [0 0.7];
    subplot(4,4,10);
    hold on
    plot(abs(sigrec{subsampling_factors(1),sig}))
    plot(abs(sigrec{1,sig}))
    ylim([0 1])

    subplot(4,4,11);
    hold on
    plot(abs(sigrec{subsampling_factors(2),sig}))
    plot(abs(sigrec{1,sig}))
    ylim([0 1])

    subplot(4,4,12);
    hold on
    plot(abs(sigrec{subsampling_factors(3),sig}))
    plot(abs(sigrec{1,sig}))
    ylim([0 1])
end

subplot(4,4,14);
hold on   
plot(interleaves, RMSE_grid(:,1,1)  *1e3, '.-', 'Color', colors(1,:))
plot(interleaves, RMSE_gllim(:,1,1) *1e3, '.-', 'Color', colors(2,:))
plot(interleaves, RMSE_nn(:,1,1) *1e3,    '.-', 'Color', colors(3,:))
plot(interleaves, RMSE_grid(:,1,2)  *1e3, 'o-', 'Color', colors(1,:))
plot(interleaves, RMSE_gllim(:,1,2) *1e3, 'o-', 'Color', colors(2,:))
plot(interleaves, RMSE_nn(:,1,2) *1e3,    'o-', 'Color', colors(3,:))
ylim([0 800])
title('T1 (ms)')

subplot(4,4,15);
hold on
plot(interleaves, RMSE_grid(:,2,1)  *1e3, '.-', 'Color', colors(1,:))
plot(interleaves, RMSE_gllim(:,2,1) *1e3, '.-', 'Color', colors(2,:))
plot(interleaves, RMSE_nn(:,2,1) *1e3,    '.-', 'Color', colors(3,:))
plot(interleaves, RMSE_grid(:,2,2)  *1e3, 'o-', 'Color', colors(1,:))
plot(interleaves, RMSE_gllim(:,2,2) *1e3, 'o-', 'Color', colors(2,:))
plot(interleaves, RMSE_nn(:,2,2) *1e3,    'o-', 'Color', colors(3,:))
ylim([0 200])
title('T2 (ms)')

if nb_param == 3
    subplot(4,4,16);
    hold on
    plot(interleaves, RMSE_grid(:,3,1)  *180/pi, '.-', 'Color', colors(1,:))
    plot(interleaves, RMSE_gllim(:,3,1) *180/pi, '.-', 'Color', colors(2,:))
    plot(interleaves, RMSE_nn(:,3,1) *180/pi,    '.-', 'Color', colors(3,:))
    plot(100,1,'k.')
    plot(100,1,'ko')
    plot(interleaves, RMSE_grid(:,3,2)  *180/pi, 'o-', 'Color', colors(1,:))
    plot(interleaves, RMSE_gllim(:,3,2) *180/pi, 'o-', 'Color', colors(2,:))
    plot(interleaves, RMSE_nn(:,3,2) *180/pi,    'o-', 'Color', colors(3,:))
    ylim([0 100])
    xlim([0 50])
    title('df (Hz)')
end
legend('DBM', 'DB-SL', 'DB-DL', ['N = ' num2str(nb_signals(1))], ['N = ' num2str(nb_signals(2))])


h(2) = subplot(442);
imagesc(1e3*T1)
colormap(h(2), hot)
colorbar
axis image; axis off

h(3) = subplot(443);
imagesc(1e3*T2)
colormap(h(3), hot)
colorbar
axis image; axis off

h(4) = subplot(444);
imagesc(df*180/pi)
colormap(h(4), hot)
colorbar
axis image; axis off


%%

if backup == 1
    savefig(fig, 'figures/AliasingNoise')
end
