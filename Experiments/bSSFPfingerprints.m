
%% Description
%
% We evaluate the DBL method with usual MRF fingerprints. This experiment
% is based on the one used to investigate the impact of dictionary size. 
%
% Fabien Boux - 01/2020

Init
disp(['Running experiment ' mfilename '.m'])


%% Setting

% Execution settings
verbose = 1;
backup  = 1;

% Signal settings
int_T1  = 1e-3 * [200 3000];	% between 20 and 3000 ms
int_T2  = 1e-3 * [20  300];     % between 20 and 3000 ms
int_df  = pi/180 * [-200 200];  % +/- 400 Hz

nb_param = 3;   % 2 for T1 and T2 estimates, and 3 for T1, T2 and Df estimates
nb_signals = [6^3 10^3];%[16^3 61^3];
S       = 50;

% Experiment settings
dbdl_computation = 1;
snr_levels = [logspace(1, 2.035, 39) inf];
nb_test_signals = 10000;

% Regression settings
Model_.K = 100;
Model_.cstr.Sigma  = 'd*';
Model_.snrtrain = 60;
%one can specify a file of previous saved models, else empty
load_model = [];


%% Creating data

% Simulation settings
st      = 5; %standard deviation for FA
FA      = (pi/180)* [90 ... %Flip anglas (rad)
           10 + 50 *sin((2*pi/500) *(1:250)) + st*randn(1,250) ...
           zeros(1,50) ...
           (10 + 50 *sin((2*pi/500) *(1:250)) + st*randn(1,250)) /2 ...
           zeros(1,50) ...
           10 + 50 *sin((2*pi/500) *(1:250)) + st*randn(1,250) ...
           zeros(1,50) ...
           (10 + 50 *sin((2*pi/500) *(1:99)) + st*randn(1,99)) /2]; % Flip angles
TR      = perlin(length(FA)+100); TR = TR(1,:); % Repetition Time (sec)
w       = gausswin(50);
TR      = filter(w,1,TR);
TR      = TR(51:1050);
TR      = 1e-3 * (10.5 + 3.5* (TR - min(TR)) ./ max(TR - min(TR)));


if ~isempty(load_model)
    load(['outputs/' load_model], 'FA','TR', 'Params', 'NeuralNet','mx','mn');
end

TR = TR(1:S);
FA = FA(1:S);


load('inputs/patterns.mat', 'FA','TR')


%% Init
t_grid          = nan(length(snr_levels), length(nb_signals));
t_gllim         = t_grid;
t_nn            = t_grid;
NRMSE_grid      = nan(length(snr_levels), length(nb_signals), nb_param);
RMSE_grid       = NRMSE_grid; MAE_grid = NRMSE_grid;
NRMSE_gllim     = NRMSE_grid;
RMSE_gllim      = NRMSE_gllim; MAE_gllim = NRMSE_gllim; 
NRMSE_nn        = NRMSE_grid;
RMSE_nn         = NRMSE_nn; MAE_nn = NRMSE_nn; 


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
    X       = (D.normalization.*D.magnetization).';
    Dico_DBM = FormatDico(X, Y(:,1:nb_param));

    
    if isempty(load_model)
        
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
        X       = (D.normalization.*D.magnetization).';
        Dico_DBL = FormatDico([real(X) imag(X)], Y(:,1:nb_param));

        if size(Dico_DBM.MRSignals,1) ~= size(Dico_DBL.MRSignals,1)
            warning('Sizes are not equals')
        end
    
        [~,Model{f}] = AnalyzeMRImages([],Dico_DBL,'DB-SL',Model_);

        if dbdl_computation == 1
            [~,NeuralNet{f}] = AnalyzeMRImages([],Dico_DBL,'DB-DL');
        end
    end
    
    clear D
    for snr = 1:length(snr_levels)

        if verbose == 2, disp(['f = ' num2str(f) ' & SNR = ' num2str(snr_levels(snr))]); end

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
        Xtest   = (D.normalization.*D.magnetization).';
        
        % Add noise
        [XtestN, tmp]       = AddNoise(Xtest, snr_levels(snr));
        real_snr(snr,f)     = mean(tmp); tmp = [];

        % Perform DBM
        Estim   = AnalyzeMRImages(XtestN,Dico_DBM,'DBM',[],Ytest(:,1:nb_param));

        t_grid(snr,f)       = Estim.GridSearch.quantification_time;
        NRMSE_grid(snr,f,:) = Estim.GridSearch.Errors.Nrmse;
        RMSE_grid(snr,f,:)	= Estim.GridSearch.Errors.Rmse;
        MAE_grid(snr,f,:) 	= Estim.GridSearch.Errors.Mae;
        
        %concat for DBL methods
        XtestN = [real(XtestN) imag(XtestN)];
        
        % Perform DB-SL
        Estim   = AnalyzeMRImages(XtestN,[],'DB-SL',Model{f},Ytest(:,1:nb_param));

        t_gllim(snr,f)   	= Estim.Regression.quantification_time;
        NRMSE_gllim(snr,f,:) = Estim.Regression.Errors.Nrmse;
        RMSE_gllim(snr,f,:) = Estim.Regression.Errors.Rmse;
        MAE_gllim(snr,f,:)  = Estim.Regression.Errors.Mae;
        
        % Perform DB-DL
        if dbdl_computation == 1
            Estim   = AnalyzeMRImages(XtestN,[],'DB-DL',NeuralNet{f},Ytest(:,1:nb_param));

            t_nn(snr,f)   	= Estim.Regression.quantification_time;
            NRMSE_nn(snr,f,:) = Estim.Regression.Errors.Nrmse;
            RMSE_nn(snr,f,:) = Estim.Regression.Errors.Rmse;
            MAE_nn(snr,f,:)  = Estim.Regression.Errors.Mae;
        end
            
    end %snr
end

t(:,1,:) 	= t_grid;
t(:,2,:) 	= t_gllim;
NRMSE(:,1,:,:) = NRMSE_grid;
NRMSE(:,2,:,:) = NRMSE_gllim;
RMSE(:,1,:,:) = RMSE_grid;
RMSE(:,2,:,:) = RMSE_gllim;
MAE(:,1,:,:) = MAE_grid;
MAE(:,2,:,:) = MAE_gllim;
if dbdl_computation == 1
    t(:,3,:) 	= t_nn;
    RMSE(:,3,:,:) = RMSE_nn;
    NRMSE(:,3,:,:) = NRMSE_nn;
    MAE(:,3,:,:) = MAE_nn;
end

fing_signals = X(1:round((1/1000)*size(X,1)):end,:);


%% Saving 

if backup == 1
    if isempty(load_model)
        save(['outputs/bSSFP_models_ID' num2str(2^16) '.mat'], 'FA','TR', 'Model', 'NeuralNet', 'nb_signals','-v7.3');
    end
    clear tmp* D Dico* X* Y*
    save(['temp/' 'bSSFPfingerprints'],'-v7.3')
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
plot(abs(fing_signals([800 30 20 1],:)') )%./ vecnorm(fing_signals([2000 1500 30 20 10],:)')) for normalized signals
plot(abs(fing_signals(randi(size(fing_signals,1),1,4),:)') )%./ vecnorm(fing_signals([2000 1500 30 20 10],:)')) for normalized signals
xlim([0 length(TR)])
xlabel('TR index');
title('(c)')


ldg = {'d','e','f'};
marks = {'.', 'o'};
param_name = {'T_1','T_2','df'};
fact = [1 1 180/pi];

for p = 1:size(RMSE,4)

    ax(p) = subplot(2,3,3+p);
    set(groot,'defaultAxesColorOrder',colors)

    hold on;
    for f = 1:length(nb_signals)
        for i = 1:size(RMSE,2)
            if f == length(nb_signals)
                plot(real_snr(:,f), squeeze(RMSE(:,i,f,p) *fact(p)), '-', 'Marker', marks{f}, 'Color', colors(i,:))
            else 
                plot(real_snr(:,f), squeeze(RMSE(:,i,f,p) *fact(p)), '-', 'Marker', marks{f}, 'Color', colors(i,:), 'HandleVisibility','off')
            end
        end
    end
    
    axis tight;
    set(gca,'FontSize',12, 'Xscale','log')
    lgd.FontSize = 16;
    
    title(param_name{p})
    
    xlabel('SNR');
    ylabel('RMSE (s)');
    
    switch p
        case 1
            ylim([0 0.6])
        case 2
            ylim([0 0.04])
        case 3
            ylim([0 120])
    end
end
plot(200, 1, 'k.')
plot(200, 1, 'ko')
legend('DBM', 'DB-SL','DB-DL', ['N = ' num2str(nb_signals(1))], ['N = ' num2str(nb_signals(2))]);


%% Exporting figures

if backup == 1
    savefig(fig, ['figures/' 'bSSFPfingerprints'])
end

