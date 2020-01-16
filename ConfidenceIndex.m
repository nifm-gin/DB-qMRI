
%% Description
%
% Explaination
%
% Fabien Boux - 01/2020


%% Setting

% Execution settings
verbose = 1;
backup  = 1;

% Signal settings
int   	= [0.1 1];
nb_param = 2;
nb_train_signals = 2000;
nb_test_signals  = 10000;
snr_train_p1 = 40;
snr_test_p1 = [20 40 60];
snr_train_p2 = [40 50 60];
snr_test_p2 = logspace(1, 2, 20);

% Experiment settings
nb_repetition_p1 = 50;
nb_repetition_p2 = 20;

% Regression settings
Parameters.K = 50;
Parameters.cstr.Sigma  = 'd*';
Parameters.cstr.Gammat = ''; 
Parameters.cstr.Gammaw = '';
Parameters.Lw = 0;


%% Creating data

% Adding to path
addpath(genpath('functions'))
addpath(genpath('tools'))

% Init
pp_std = nan(length(snr_test_p1),nb_param,nb_repetition_p1);
pp_err = pp_std;


%% Processing 1

disp('First process')

parfor rep = 1:nb_repetition_p1
    
    if verbose == 1, disp([num2str(rep) '/' num2str(nb_repetition_p1)]); end
    
    p   = [.01 .01];
    while min(abs(pdist(p',@(x,y) x-y))) < .1
        p   = 0.1 + 0.9*rand(1,5);
    end

    % Generate training signals
    Xtrain  = [];
    Ytrain	= int(1) + (int(2)-int(1)) * net(scramble(sobolset(nb_param),'MatousekAffineOwen'),nb_train_signals);
    for sim = 1:size(Ytrain,1)
        Xtrain(sim,:) = toyMRsignal(Ytrain(sim,:),p(1:nb_param));
    end

    % Generate test signals
    Xtest   = [];
    Ytest	= int(1) + (int(2)-int(1)) * rand(nb_test_signals,nb_param);
    for sim = 1:size(Ytest,1)
        Xtest(sim,:) = toyMRsignal(Ytest(sim,:),p(1:nb_param));
    end

    % GlliM learning
    Dico = [];
    Dico{1}.MRSignals       = AddNoise(Xtrain, snr_train_p1);
    Dico{1}.Parameters.Par  = Ytrain;
    [~,Params]  = AnalyzeMRImages(Xtrain,Dico,'RegressionMRF',Parameters);
    
    % Estimates
    tmp_pp_std = [];
    tmp_pp_err = [];
    for s = 1:length(snr_test_p1)
        
        Xtest_noisy = AddNoise(Xtest,snr_test_p1(s));
        
        Estim   = AnalyzeMRImages(Xtest_noisy, [], 'RegressionMRF', Params);
        Ygllim  = squeeze(Estim.Regression.Y(:,1:nb_param));
        Cov     = squeeze(Estim.Regression.Cov);
        Rmse    = EvaluateEstimation(Ytest, Ygllim);

        tmp_pp_std(s,:) = nanmean(Cov,1).^.5;
        tmp_pp_err(s,:) = nanmean((Ytest-Ygllim).^2,1).^.5;
    end
    
    pp_std(:,:,rep) = tmp_pp_std;
    pp_err(:,:,rep) = tmp_pp_err;
end

%% Saving (process 1)

if backup == 1
    clear tmp* X* Y* Dico* Estim
    save(['temp/' 'ConfidenceIndex-1'])
end


%% Processing 2

disp('Second process')

cc  = nan(nb_param, length(snr_test_p2), nb_repetition_p2, length(snr_train_p2));

for st = 1:length(snr_train_p2)
    
    tmp_pp_std = nan(nb_param, length(snr_test_p2));
    tmp_pp_err = tmp_pp_std;

    if verbose == 1, disp(['SNR_train ' num2str(st) '/' num2str(length(snr_train_p2))]); end
    
    parfor rep = 1:nb_repetition_p2

        if verbose > 1, disp([num2str(rep) '/' num2str(nb_repetition_p2)]); end

        p   = [.01 .01];
        while min(abs(pdist(p',@(x,y) x-y))) < .1
            p   = 0.1 + 0.9*rand(1,5);
        end

        % Generate training signals
        Xtrain  = [];
        Ytrain	= int(1) + (int(2)-int(1)) * net(scramble(sobolset(nb_param),'MatousekAffineOwen'),nb_train_signals);
        for sim = 1:size(Ytrain,1)
            Xtrain(sim,:) = toyMRsignal(Ytrain(sim,:),p(1:nb_param));
        end
        [Xtrain,snrtmp] = AddNoise(Xtrain, snr_train_p2(st));
        real_snr_train(:,rep) = nanmean(snrtmp);

        % Generate test signals
        Xtest   = [];
        Ytest	= int(1) + (int(2)-int(1)) * rand(nb_test_signals,nb_param);
        for sim = 1:size(Ytest,1)
            Xtest(sim,:) = toyMRsignal(Ytest(sim,:),p(1:nb_param));
        end

        % GlliM learning
        Dico    = [];
        Dico{1}.MRSignals       = Xtrain;
        Dico{1}.Parameters.Par  = Ytrain;
        [~,Params]  = AnalyzeMRImages(Xtrain,Dico,'RegressionMRF',Parameters);

        % Estimating
        tmp_pp_std = [];
        tmp_pp_err = [];
        for s = 1:length(snr_test_p2)

            Xtest_noisy = AddNoise(Xtest, snr_test_p2(s));

            Estim   = AnalyzeMRImages(Xtest_noisy, [], 'RegressionMRF', Params);
            Ygllim  = squeeze(Estim.Regression.Y(:,1:nb_param));
            Cov     = squeeze(Estim.Regression.Cov);
            Rmse    = EvaluateEstimation(Ytest, Ygllim);

            tmp_pp_std(:,s) = nanmean(Cov,1).^.5;
            tmp_pp_err(:,s) = nanmean((Ytest-Ygllim).^2,1).^.5;
        end

        cc(:,:,rep,st) = tmp_pp_std ./ tmp_pp_err;
    end
end


%% Saving (process 2)

if backup == 1
    clear tmp* X* Y* Dico* Estim
    save(['temp/' 'ConfidenceIndex-2'])
end


%% Displaying

fig = figure;

subplot(121)
hold on

for f = 1:size(pp_std,1)
    
    plot(squeeze(mean(pp_std(f,:,:),2))', squeeze(mean(pp_err(f,:,:),2))','.','MarkerSize',12)

    mdl = fitlm(squeeze(mean(pp_std(f,:,:),2)), squeeze(mean(pp_err(f,:,:),2)),'Intercept',false);
    plot([0 max(pp_std(:))],[0 max(pp_std(:))]*mdl.Coefficients.Estimate, 'k')

    leg{2*f-1}  = ['SNR_{test} = ' num2str(snr_test_p1(f))];
    leg{2*f}    = ['\alpha = ' num2str(mdl.Coefficients.Estimate,3) '   (R^2 = ' num2str(mdl.Rsquared.Ordinary,2) ')'];
end

title('RMSE = \alpha x CI')
xlabel('CI'); ylabel('RMSE')
legend(leg)


subplot(122)
hold on
clear leg

F = @(xdata,x)x(1)*exp(-x(2)./xdata) + x(3);

for st = 1:length(snr_train_p2)
    
    ff = mean(mean(1./squeeze(cc(:,:,:,st)),1),2);
    plot(snr_test_p2, ff,'.-','Markersize',15)

    a = levenbergmarquardt(F, snr_test_p2, ff', [max(ff) -1 min(ff)]);
    plot(snr_test_p2, F(snr_test_p2,a), 'k--')

    leg{2*st-1} = ['SNR_{train} = ' num2str(snr_train_p2(st))];
    leg{2*st} = 'fit';
end
plot([10 100], [1 1], 'k')

title('\alpha = f(SNR_{test})')
xlabel('SNR_{test}'); ylabel('\alpha')
legend(leg)


%% Exporting figures

if backup == 1
    savefig(fig, ['figures/' 'ConfidenceIndex'])
end



