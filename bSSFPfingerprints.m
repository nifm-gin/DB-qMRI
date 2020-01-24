
%% Description
%
% We evaluate the DBL method with usual MRF fingerprints. This experiment
% is based on the one used to investigate the impact of dictionary size. 
%
% Fabien Boux - 01/2020


%% Setting

% Execution settings
verbose = 1;
backup  = 1;

% Signal settings
int_T1  = 1e-3 * [100 5100];	% between 100 and 5000 ms
int_T2  = 1e-3 * [20  3020];    % between 20 and 3000 ms
int_df  = pi/180 * [-400 400];  	% +/- 400 Hz

nb_param = 3;   % 2 for T1 and T2 estimates, and 3 for T1, T2 and Df estimates
% nb_signals = [400 6400]; % if nb_param == 2
nb_signals = [4096 132651]; % if nb_param == 3
    
% Simulation settings
FA      = (pi/180)* [90 ...
           10 + 50 *sin((2*pi/500) *(1:250)) + 5*randn(1,250) ...
           zeros(1,49) ...
           (10 + 50 *sin((2*pi/500) *(1:200)) + 5*randn(1,200)) /2]; % Flip angles
% TR      = 1e-3 * (10.5  + 3.5 * rand(1,500));    % Uniform between 10.5 and 14 ms

% Experiment settings
snr_levels = [logspace(1, 2.035, 39) inf];
nb_test_signals = 10000;

% Regression settings
Parameters.K = 100;
Parameters.cstr.Sigma  = 'd*';
Parameters.cstr.Gammat = ''; 
Parameters.cstr.Gammaw = '';
Parameters.Lw = 0;
snr_train = 60;


%% Creating data

% Adding to path
addpath(genpath('functions'))
addpath(genpath('tools'))

%
TR      = perlin(500); TR = TR(1,:);
TR      = (TR - min(TR)) ./ max(TR - min(TR));
TR      = 1e-3 * (10.5  + 3.5 * TR);

% Init
t_grid          = nan(length(snr_levels), length(nb_signals));
t_gllim         = t_grid;
mNRMSE_grid     = nan(length(snr_levels), length(nb_signals), nb_param);
mRMSE_grid      = mNRMSE_grid; mMAE_grid = mNRMSE_grid;
mNRMSE_gllim    = mNRMSE_grid;
mRMSE_gllim     = mNRMSE_gllim; mMAE_gllim = mNRMSE_gllim; 


%% Processing 

if verbose >= 1, disp('Dictionary size'); end

for f = 1:size(nb_signals,2)

    if verbose >= 1, disp(num2str(nb_signals(f))); end

    % Compute dico grid
    clear X Y
    nb_step = round(nb_signals(f)^(1/nb_param));
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
    X       = abs(D.normalization.*D.magnetization)';
    DicoG{1}.MRSignals = X; 
    DicoG{1}.Parameters.Par = Y(:,1:nb_param);

    % Compute training dataset
    clear X Y
    Y       = net(scramble(sobolset(nb_param),'MatousekAffineOwen'),nb_signals(f));
    Y(:,1)  = int_T1(1) + (int_T1(2) - int_T1(1)) * Y(:,1);
    Y(:,2)  = int_T2(1) + (int_T2(2) - int_T2(1)) * Y(:,2); 
    if nb_param == 2
        D       = MRF_dictionary(Y(:,1), Y(:,2), [], FA, TR); 
    elseif nb_param == 3
        Y(:,3)  = int_df(1) + (int_df(2) - int_df(1)) * Y(:,3);
        D       = MRF_dictionary(Y(:,1), Y(:,2), Y(:,3), FA, TR); 
    end 
    X       = abs(D.normalization.*D.magnetization)';
    DicoR{1}.MRSignals = AddNoise(X, snr_train); 
    DicoR{1}.Parameters.Par = Y(:,1:nb_param);

    if size(DicoG{1}.MRSignals,1) ~= size(DicoR{1}.MRSignals,1)
        warning('Sizes are not equals')
    end
    
    [~, Params] = AnalyzeMRImages([],DicoR,'DBL',Parameters);
    
    clear D
    parfor snr = 1:length(snr_levels)

        if verbose == 2, disp(['f = ' num2str(f) ' Snr order: ' num2str(snr_levels(snr))]); end

        % Generate test data
        Ytest 	= rand(nb_test_signals,nb_param);
        Ytest(:,1) = int_T1(1) + (int_T1(2) - int_T1(1)) * Ytest(:,1);
        Ytest(:,2) = int_T2(1) + (int_T2(2) - int_T2(1)) * Ytest(:,2); 
        if nb_param == 2
            D       = MRF_dictionary(Ytest(:,1), Ytest(:,2), [], FA, TR); 
        elseif nb_param == 3
            Ytest(:,3) = int_df(1) + (int_df(2) - int_df(1)) * Ytest(:,3);
            D       = MRF_dictionary(Ytest(:,1), Ytest(:,2), Ytest(:,3), FA, TR); 
        end
        Xtest   = abs(D.normalization.*D.magnetization)';
        
        % Add noise
        [XtestN, tmp]       = AddNoise(Xtest, snr_levels(snr));
        real_snr(snr,f)     = mean(tmp); tmp = [];

        % Perform DBM
        tic;
        Estim   = AnalyzeMRImages(XtestN,DicoG,'DBM',[],Ytest(:,1:nb_param));

        t_grid(snr,f)       = toc;
        mNRMSE_grid(snr,f,:) = Estim.GridSearch.Errors.Nrmse;
        mRMSE_grid(snr,f,:)	= Estim.GridSearch.Errors.Rmse;
        mMAE_grid(snr,f,:) 	= Estim.GridSearch.Errors.Mae;

        % Perform DBL
        tic;
        Estim   = AnalyzeMRImages(XtestN,[],'DBL',Params,Ytest(:,1:nb_param));

        t_gllim(snr,f)   	= toc;
        mNRMSE_gllim(snr,f,:) = Estim.Regression.Errors.Nrmse;
        mRMSE_gllim(snr,f,:) = Estim.Regression.Errors.Rmse;
        mMAE_gllim(snr,f,:) = Estim.Regression.Errors.Mae;

    end %snr
end

t(:,1,:) 	= t_grid;
t(:,2,:) 	= t_gllim;
mNRMSE(:,1,:,:) = mNRMSE_grid;
mNRMSE(:,2,:,:) = mNRMSE_gllim;
mRMSE(:,1,:,:) = mRMSE_grid;
mRMSE(:,2,:,:) = mRMSE_gllim;
mMAE(:,1,:,:) = mMAE_grid;
mMAE(:,2,:,:) = mMAE_gllim;

fing_signals    = X(randi(size(X,1),10,1),:);


%% Saving 

if backup == 1
    clear tmp* D Dico* X* Y*
    save(['temp/' 'bSSFPfingerprints'])
end


%% Displaying

colors  = [          0    0.4470    0.7410
                0.8500    0.3250    0.0980
                0.9290    0.6940    0.1250
                0.4940    0.1840    0.5560
                0.4660    0.6740    0.1880
                0.3010    0.7450    0.9330
                0.6350    0.0780    0.1840];
colors = [colors; colors];
    
fig = figure;

subplot(4,4,[1 2])
plot(FA * 180 / pi)
xlabel('TR index'); ylabel('FA (degrees)')
ylim([0 90]); xlim([0 length(TR)])
title('(a)')

subplot(4,4,[5 6])
plot(TR*1e3)
xlabel('TR index'); ylabel('TR (ms)')
ylim([10 14]); xlim([0 length(TR)])
title('(b)')

subplot(4,4,[3 4 7 8])
plot(fing_signals')
% ylim([0 .2]); 
xlim([0 length(TR)])
xlabel('TR index');
title('(c)')


ldg = {'d','e','f'};

count = 1;
for f = 1:length(nb_signals)

    ax(f) = subplot(2,2,2+f);
    set(groot,'defaultAxesColorOrder',colors)

    hold on;
%     plot([0 200], [mRMSE(end,1,f,1) mRMSE(end,1,f,1)], '--', 'linewidth',1.5, 'Color',colors(1,:),'HandleVisibility','off')        
%     plot([0 200], [mRMSE(end,1,f,2) mRMSE(end,1,f,2)], '--', 'linewidth',1.5, 'Color',colors(2,:),'HandleVisibility','off')        

    plot(real_snr(:,f), squeeze(mRMSE(:,1,f,1)), '*-', 'MarkerSize', 12, 'Color', colors(1,:))
    plot(real_snr(:,f), squeeze(mRMSE(:,1,f,2)), '*-', 'MarkerSize', 12, 'Color', colors(2,:))

%     plot([0 200], [mRMSE(end,2,f,1) mRMSE(end,2,f,1)], '--', 'linewidth',1.5, 'Color',colors(1,:),'HandleVisibility','off')        
%     plot([0 200], [mRMSE(end,2,f,2) mRMSE(end,2,f,2)], '--', 'linewidth',1.5, 'Color',colors(2,:),'HandleVisibility','off')        

    plot(real_snr(:,f), squeeze(mRMSE(:,2,f,1)), '.-', 'MarkerSize', 12, 'Color', colors(1,:))
    plot(real_snr(:,f), squeeze(mRMSE(:,2,f,2)), '.-', 'MarkerSize', 12, 'Color', colors(2,:))
    
    if f ==  size(nb_signals,2), legend('T_1 - DBM','T_2 - DBM', 'T_1 - DBL','T_2 - DBL'); end
    title([num2str(nb_signals(f)) ' signals']);
    if count == 1, ylabel('RMSE (s)'); end

    axis tight;
    set(gca,'FontSize',12)
    lgd.FontSize = 16;

    title(['(' ldg{f} ') ' num2str(nb_signals(f)) '-signal dictionary'])
    
%     xlim([5 110])
%     ylim([0.01 0.5])

    xlabel('SNR');
end
linkaxes(ax,'y')


%% Exporting figures

if backup == 1
    savefig(fig, ['figures/' 'bSSFPfingerprints-supp'])
end

