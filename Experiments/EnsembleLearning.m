
%% Description
%
% We investigate a combination of dictionary-based learning methods (ie 
% DB-SL and DB-DL) for robust estimation
%
% Fabien Boux - 11/2020

Init
disp(['Running experiment ' mfilename '.m'])


%%

exec_nn = 'cpu';

% Signal settings
int     = [.01 1];
p       = [.01 .01];
while min(abs(pdist(p',@(x,y) x-y))) < .05, p = 0.1 + 0.9*rand(1,10); end
p = [0.9996    0.6002    0.2745    0.8402    0.4683    0.7702    0.6516    0.1095    0.3859    0.9407];

% Settings
nb_param = 3;
nb_signals = 4e3;
nb_test_signals = 1e4;
%snr = [20 30 40 80];
snr = [logspace(1, 2.035, 39) inf];

% Regression settings
Parameters = [];
Parameters.K = 50;
Parameters.cstr.Sigma  = 'd*';
Parameters.cstr.Gammat = ''; 
Parameters.cstr.Gammaw = '';
Parameters.Lw = 0;
snr_train = 60;


% Compute training dataset
Ytrain  	= int(1) + (int(2) - int(1)) * net(scramble(sobolset(nb_param),'MatousekAffineOwen'), nb_signals);
parfor sim = 1:size(Ytrain,1)
    Xtrain(sim,:) = toyMRsignal(Ytrain(sim,:),p(1:nb_param));
end
DicoR{1}.MRSignals = AddNoise(abs(Xtrain), snr_train); 
DicoR{1}.Parameters.Par = Ytrain;


Ytest 	= int(1) + (int(2) - int(1)) * rand(nb_test_signals, nb_param);
Xtest = [];
parfor sim = 1:size(Ytest,1)
    Xtest(sim,:) = toyMRsignal(Ytest(sim,:),p(1:nb_param));
end

[Estim,Params] = AnalyzeMRImages(DicoR{1}.MRSignals, DicoR, 'DBL', Parameters);
Y = squeeze(Estim.Regression.Y);
CI_ = squeeze(Estim.Regression.Cov.^.5);
disp('GLLiM learned')

mn      = min(DicoR{1}.Parameters.Par);
mx      = max(DicoR{1}.Parameters.Par - mn);
YtrainNN = (DicoR{1}.Parameters.Par - mn) ./ mx;
NeuralNet = EstimateNNmodel(DicoR{1}.MRSignals,YtrainNN,0,exec_nn);
disp('NN learned')

%
% NeuralNet_comb_ci = EstimateNNmodel([Ynn_ Y CI_], YtrainNN, 0, exec_nn, 50,6);
NeuralNet_comb_ci = EstimateNNmodel([Y CI_], YtrainNN, 0, exec_nn, 50,6);
disp('Combined NN with CI learned')

%
% NeuralNet_comb = EstimateNNmodel([Ynn_ Y], YtrainNN, 0, exec_nn, 50,6);
NeuralNet_comb = EstimateNNmodel([Y], YtrainNN, 0, exec_nn, 50,6);
disp('Combined NN learned')


%%

for i = 1:numel(snr)
    
    [XtestN, tmp] = AddNoise(Xtest, snr(i));
    real_snr(i) = mean(tmp);
    
    %gllim
    Estim   = AnalyzeMRImages(XtestN, [], 'DBL', Params, [],[], snr(i));
    Ygllim  = squeeze(Estim.Regression.Y);
    CI      = squeeze(Estim.Regression.Cov.^.5);

    %nn
    Ynn     = EstimateParametersFromNNmodel(XtestN,NeuralNet,exec_nn);
    Ynn    	= (Ynn .* mx) + mn;
    
    %mean
    Ymean = (Ygllim + Ynn) / 2;
    
    %mean ci
    thresh = max(CI(:));
    p = 1 - CI / thresh;
    p(p < 0) = 0;
    Yci = (Ygllim .* p) + (Ynn .* (1-p));
    
    
    %comb
    Ycomb 	= EstimateParametersFromNNmodel([Ygllim],NeuralNet_comb,exec_nn);
    Ycomb  	= (Ycomb .* mx) + mn;
    
    %comb ci
    Ycomb_ci = EstimateParametersFromNNmodel([Ygllim CI],NeuralNet_comb_ci,exec_nn);
    Ycomb_ci = (Ycomb_ci .* mx) + mn;
    
    
    Rmse_gllim(i) = mean(EvaluateEstimation(Ytest, Ygllim));
    Rmse_nn(i) = mean(EvaluateEstimation(Ytest, Ynn));
    Rmse_mean(i) = mean(EvaluateEstimation(Ytest, Ymean));
    Rmse_mean_ci(i) = mean(EvaluateEstimation(Ytest, Yci));
    Rmse_comb(i) = mean(EvaluateEstimation(Ytest, Ycomb));
    Rmse_comb_ci(i) = mean(EvaluateEstimation(Ytest, Ycomb_ci));
end


%%

fig = figure;
h(1) = subplot(121);
hold on
plot(real_snr, Rmse_gllim)
plot(real_snr, Rmse_nn)
plot(real_snr, Rmse_mean)
plot(real_snr, Rmse_mean_ci)
legend('GLLiM','NN','Mean','Mean CI')
ylabel('RMSE (s)')
xlabel('SNR')

set(gca, 'Xscale','log','Yscale','log')

h(2) = subplot(122);
hold on
plot(real_snr, Rmse_gllim)
plot(real_snr, Rmse_nn)
plot(real_snr, Rmse_comb)
plot(real_snr, Rmse_comb_ci)
legend('GLLiM', 'NN', 'Comb','Comb CI')

set(gca, 'Xscale','log','Yscale','log')






