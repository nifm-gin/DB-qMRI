
%% Description
%
% We investigate the correlation between the confidence index (CI) and the
% estimate errors (RMSE).
%
% Fabien Boux - 02/2020


%% Setting

% Execution settings
verbose = 1;
backup  = 1;

% Signal settings
int   	= [0.01 1];
nb_param = 3;
nb_train_signals = 10000;
nb_test_signals  = 10000;
snr_train = 60;
snr_test = [20 30 40 60 100];

% Experiment settings
nb_repetition = 100;

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
pp_std = nan(length(snr_test),nb_param,nb_repetition);
pp_err = pp_std; pp_err2 = pp_std; pp_std2 = pp_err2;


%% Processing

parfor rep = 1:nb_repetition
    
    if verbose == 1, disp([num2str(rep) '/' num2str(nb_repetition)]); end
    
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
    [Dico{1}.MRSignals, real_snr_train] = AddNoise(Xtrain, snr_train);
    Dico{1}.Parameters.Par  = Ytrain;
    [~,Params]  = AnalyzeMRImages(Xtrain,Dico,'DBL',Parameters);
    
    % Estimates
    tmp_pp_std = [];
    tmp_pp_err = [];
    tmp_pp_std2 = [];
    tmp_pp_err2 = [];
    for s = 1:length(snr_test)
        
        [Xtest_noisy, real_snr] = AddNoise(Xtest, snr_test(s));
        
        Estim   = AnalyzeMRImages(Xtest_noisy, [], 'DBL', Params);
        Ygllim  = squeeze(Estim.Regression.Y(:,1:nb_param));
        Cov     = squeeze(Estim.Regression.Cov);
        Rmse    = EvaluateEstimation(Ytest, Ygllim);

        tmp_pp_std(s,:) = nanmean(Cov,1).^.5;
        tmp_pp_err(s,:) = Rmse;
        
        Params_updt = Params;
        var_noise = mean((max(Xtest,[],2) ./ real_snr).^2);
        Params_updt.theta = updateSigma(Params.theta, var_noise);
        Estim   = AnalyzeMRImages(Xtest_noisy, [], 'DBL', Params_updt);
        
        %this line is equivalent to the previous implementation since it
        %corresponds to the integration in our estimation function of the
        %model correction.
%         Estim   = AnalyzeMRImages(Xtest_noisy, [], 'DBL', Params, [],[], real_snr);

        Ygllim  = squeeze(Estim.Regression.Y(:,1:nb_param));
        Cov     = squeeze(Estim.Regression.Cov);
        Rmse    = EvaluateEstimation(Ytest, Ygllim);

        tmp_pp_std2(s,:) = nanmean(Cov,1).^.5;
        tmp_pp_err2(s,:) = Rmse;
    end
    
    pp_std(:,:,rep) = tmp_pp_std;
    pp_err(:,:,rep) = tmp_pp_err;
    
    pp_std2(:,:,rep) = tmp_pp_std2;
    pp_err2(:,:,rep) = tmp_pp_err2;
end


%% Saving

if backup == 1
    clear tmp* X* Y* Dico* Estim
    save(['temp/' 'ConfidenceIndex'])
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

param = 1; % first or second parameter lead to ~ same results

subplot(131)
hold on

for f = 1:size(pp_std,1)
    plot(squeeze(pp_std(f,param,:))', squeeze(pp_err(f,param,:))', '.','MarkerSize',12, 'color',colors(f,:))

    mdl = fitlm(squeeze(pp_std(f,param,:)), squeeze(pp_err(f,param,:)),'Intercept',false);
    plot([-1 2*max(pp_std(:))], [-1 2*max(pp_std(:))]*mdl.Coefficients.Estimate, 'k')

    leg{2*f-1}  = ['SNR_{test} = ' num2str(snr_test(f))];
    leg{2*f}    = ['\alpha = ' num2str(mdl.Coefficients.Estimate,3) '   (R^2 = ' num2str(mdl.Rsquared.Ordinary,3) ')'];
end
xlim([0 max(pp_std2(:))]); ylim([0 max(pp_std2(:))])
title('RMSE = \alpha x CI')
xlabel('CI'); ylabel('RMSE')
legend(leg)


clear leg
subplot(132)
hold on
for f = 1:size(pp_std2,1)
    plot(squeeze(pp_std2(f,param,:))', squeeze(pp_err2(f,param,:))', '.','MarkerSize',12, 'color',colors(f,:))

    leg{f}  = ['SNR_{test} = ' num2str(snr_test(f))];
end
mdl = fitlm(reshape(pp_std2(:,param,:),1,[]), reshape(pp_err2(:,param,:),1,[]),'Intercept',false);
plot([-1 2*max(pp_std2(:))], [-1 2*max(pp_std2(:))]*mdl.Coefficients.Estimate, 'k--')
xlim([0 max(pp_std2(:))]); ylim([0 max(pp_std2(:))])
leg{f+1} = ['\alpha = ' num2str(mdl.Coefficients.Estimate,3) '   (R^2 = ' num2str(mdl.Rsquared.Ordinary,3) ')'];

title('RMSE_{corrected} = \alpha x CI_{corrected}')
xlabel('CI corrected'); ylabel('RMSE with Sigma corrected')
legend(leg)


clear leg
subplot(133)
hold on
for f = 1:size(pp_std2,1)
    plot(squeeze(pp_err(f,param,:))', squeeze(pp_err2(f,param,:))', '.','MarkerSize',12, 'color',colors(f,:))
    
    leg{f}  = ['SNR_{test} = ' num2str(snr_test(f))];
end
plot([-1 2*max(pp_err(:))], [-1 2*max(pp_err(:))], 'k')
xlim([0 max(pp_std2(:))]); ylim([0 max(pp_std2(:))])
leg{f+1} = 'RMSE_{corrected} = RMSE';

title('RMSE_{corrected} vs. RMSE')
xlabel('RMSE'); ylabel('RMSE with Sigma corrected')
legend(leg)


%% Exporting figures

if backup == 1
    savefig(fig, ['figures/' 'ConfidenceIndex'])
end
