
p = [1 2];

K = 100;

snr_levels = 100;
snr_train = 60;


%% Wang sequence

load('/home_eq/bouxfa/Code/MRVox2D/outputs/test/Dico_WangSeq.mat')

r = randperm(size(Dico.MRSignals,1),size(Dico.MRSignals,1));
Dico.MRSignals = Dico.MRSignals(r,:);
Dico.Parameters.Par = Dico.Parameters.Par(r,p);

samples = 1:1000;

Xtrain = Dico.MRSignals(1:4000,samples);
Ytrain = Dico.Parameters.Par(1:4000,:);
Xtest = Dico.MRSignals(4001:end,samples);
Ytest = Dico.Parameters.Par(4001:end,:);
clear Dico

%remove nan
v = any(isnan(Xtrain)') | any(isnan(Ytrain)');
Xtrain(v,:) = [];
Ytrain(v,:) = [];
v = any(isnan(Xtest)') | any(isnan(Ytest)');
Xtest(v,:) = [];
Ytest(v,:) = [];

Dico{1}.MRSignals = Xtrain;
Dico{1}.Parameters.Par = Ytrain;

% Model tuning
Parameters = [];
Parameters.K = K;
Parameters.cstr.Sigma  = 'd*';
Parameters.cstr.Gammat = ''; 
Parameters.cstr.Gammaw = '';
Parameters.Lw = 0;

% Add noise
XtestN = AddNoise(Xtest, snr_levels);

% Perform DBM
Estim = AnalyzeMRImages(XtestN,Dico,'DBM',[],Ytest);
RMSE_grid = Estim.GridSearch.Errors.Rmse;

% Perform DBL
Estim = AnalyzeMRImages(XtestN,Dico,'DBL',Parameters,Ytest);
RMSE_gllim = Estim.Regression.Errors.Rmse;


%% Wang sequence

load('/home_eq/bouxfa/Code/MRVox2D/outputs/test/Dico_MGEFIDSE.mat')

r = randperm(size(Dico.MRSignals,1),size(Dico.MRSignals,1));
Dico.MRSignals = Dico.MRSignals(r,:);
Dico.Parameters.Par = Dico.Parameters.Par(r,p);

Xtrain = Dico.MRSignals(1:4000,:);
Ytrain = Dico.Parameters.Par(1:4000,:);
Xtest = Dico.MRSignals(4001:end,:);
Ytest = Dico.Parameters.Par(4001:end,:);
clear Dico

%remove nan
v = any(isnan(Xtrain)') | any(isnan(Ytrain)');
Xtrain(v,:) = [];
Ytrain(v,:) = [];
v = any(isnan(Xtest)') | any(isnan(Ytest)');
Xtest(v,:) = [];
Ytest(v,:) = [];

Dico{1}.MRSignals = Xtrain;
Dico{1}.Parameters.Par = Ytrain;

% Model tuning
Parameters = [];
Parameters.K = K;
Parameters.cstr.Sigma  = 'd*';
Parameters.cstr.Gammat = ''; 
Parameters.cstr.Gammaw = '';
Parameters.Lw = 0;

% Add noise
XtestN = AddNoise(Xtest, snr_levels);

% Perform DBM
Estim = AnalyzeMRImages(XtestN,Dico,'DBM',[],Ytest);
RMSE_grid = Estim.GridSearch.Errors.Rmse;

% Perform DBL
Estim = AnalyzeMRImages(XtestN,Dico,'DBL',Parameters,Ytest);
RMSE_gllim = Estim.Regression.Errors.Rmse;


%%

disp('Sequence: Wang / MGEFIDSE')

for i = 1:2
    disp(' ')
    disp(['Parameter ' num2str(i)])
    disp(['DBM: ' num2str(RMSE_grid(i),3) ' / ' num2str(RMSE_grid_vmrf(i),3)]) 
    disp(['DBL: ' num2str(RMSE_gllim(i),3) ' / ' num2str(RMSE_gllim_vmrf(i),3)]) 
end