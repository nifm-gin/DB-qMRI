
%% Description
%
% TODO
%
% Fabien Boux - 11/2020

Init
disp(['Running experiment ' mfilename '.m'])


%% Setting

% Execution settings
verbose = 1;
backup  = 1;

% Experiment settings
nb_test_signals = 5e4;
nb_signals = [2e3 1e4];

K   = 3:3:200;
constraints = {'d*', 'i*'};

% Regression settings
Model.Lw = 0;
Parameters.cstr.Gammat = ''; 
Parameters.cstr.Gammaw = '';
Model.maxiter = 250;
Model.snrtrain = 60;


%% Load data set

%dictionary filenames
dicos   = {'inputs/3DvMRF/qMC.mat'}; %qMC
      
% Load echotimes
echtime = load('inputs/echoetime.mat');
echtime = echtime.t;

% Quasi-random
load(dicos{1})
Xqmc    = Dico.MRSignals(:,end/2+1:end) ./ Dico.MRSignals(:,1:end/2);
Xqmc 	= interp1(Dico.Tacq(1:size(Xqmc,2))', Xqmc', echtime)';
Yqmc    = Dico.Parameters.Par;

Xqmc(Dico.Parameters.Par(:,4)== 1e-6,:) = [];
Yqmc(Dico.Parameters.Par(:,4)== 1e-6,:) = [];

Yqmc  = Yqmc(:, [2 3 1]);

[rows, ~] = find(~isnan(Xqmc));
r       = unique([rows; find(any(isnan(Yqmc),2) == 0)]);
Xqmc    = Xqmc(r,:);
Yqmc    = Yqmc(r,:);

Xtest   = Xqmc(nb_signals(end)+1:end,:);
Ytest   = Yqmc(nb_signals(end)+1:end,:);

% Add noise
Xqmc    = AddNoise(Xqmc, Model.snrtrain);
Xtest   = AddNoise(Xtest, 50);
        

%% Creating data

% Err = struct(length(nb_params),length(constraints));
% Bic = struct(length(nb_params),length(constraints));

for c1 = 1:length(constraints)
    for c2 = 1:length(nb_signals)

        % Choose number of parameters and constrint on covariance matrices
        Model.cstr.Sigma = constraints{c1};

        % Generate training dataset
        Xtrain  = Xqmc(1:nb_signals(c2),:);
        Ytrain  = Yqmc(1:nb_signals(c2),:);

        % Init
        Rmse    = nan(length(K), size(Ytrain,2));
        bic     = nan(1,length(K));
        [N,Lt]  = size(Ytrain);
       	Lw      = Model.Lw;
        D       = size(Xtrain,2);

        % Start evaluating K impact
        parfor k = 1:length(K)
            disp(k)
            
            try
                v = randi(size(Xtest,1),nb_test_signals,1);
                
                [theta,~,ll] = EstimateInverseFunction(Ytrain, Xtrain, K(k), Model.Lw, Model.maxiter, Model.cstr, 0);
                Ypredict	= EstimateParametersFromModel(Xtest(v,:), theta, 0);

                Rmse(k,:)   = EvaluateEstimation(Ytest(v,:), Ypredict);

                switch Model.cstr.Sigma
                    case 'i*'
                        M = K(k)* (D*(Lw+Lt+1)     + Lt*(Lt+3)/2 + 1);
                    case 'i'
                        M = K(k)* (D*(Lw+Lt+1)     + Lt*(Lt+3)/2 + 1)  + K(k)-1;
                    case 'd*'
                        M = K(k)* (D*(Lw+Lt+2)     + Lt*(Lt+3)/2) + D  + K(k)-1;
                    case 'd'
                        M = K(k)* (D*(Lw+Lt+2)     + Lt*(Lt+3)/2)      + K(k)-1;
                    otherwise
                        M = K(k)* (D*(Lw+Lt+D+1)   + Lt*(Lt+3)/2)      + K(k)-1;
                end

                bic(k) = -2*ll + M*log(N);
            
            catch
                warning(['Iteration ' Model.cstr.Sigma '/' num2str(K(k)) ' failed'])
            end
        end

        Err{c2,c1}  = Rmse;
        Bic{c2,c1}  = bic;
    end
end


%% Backup

if backup == 1
    clear tmp* Dico* Estim X* Y*
    save(['temp/' 'Kimpact_vMRF'])
end


%% Disp

fig = figure;
hold on
for p = 1:3
    subplot(1,3,p)
    hold on
    
    for c2 = 1:length(nb_signals)
        plot(K, mean(Err{c2,1}(:,p),2), '.-', 'linewidth',1.5);
        %plot(K, mean(Dic{c2,c1},2), '.-', 'linewidth',1.5);
    end
end
ylabel('RMSE')
%ylabel('BIC')
xlabel('K')
% legend([repelem('P = ',length(nb_params),1) num2str(nb_params')])


%%

if backup == 1
    savefig(fig, ['figures/' 'Kimpact_vMRF'])
end