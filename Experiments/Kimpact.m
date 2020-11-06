
%% Description
%
% The DBL method estimates the parameters using a continuous function that
% is not restricted to the dictionary parameter space. We therefore 
% investigated the behaviors of the DBL and DBM methods outside the
% boundaries of this space. Two 2-dimensional parameter subspaces was
% defined: one to sample the parameters of the dictionary, and a larger one
% for the test.
%
% Fabien Boux - 01/2020

Init
disp(['Running experiment ' mfilename '.m'])


%% Setting

% Execution settings
verbose = 1;
backup  = 1;

% Signal settings
int     = [.01 1];
p       = [.01 .01];
while min(abs(pdist(p',@(x,y) x-y))) < .05, p = 0.1 + 0.9*rand(1,10); end

% Experiment settings
nb_test_signals = 1e4;
nb_train_signals = 1e4;

K = 3:3:200;
constraints = {'d*'}; %{'i*','i','d*','d',''};
nb_params = [3 4 5 7];

% Regression settings
Model.Lw = 0;
Parameters.cstr.Gammat = ''; 
Parameters.cstr.Gammaw = '';
Model.maxiter = 250;
Model.snrtrain = 60;


%% Creating data

Err = struct(length(nb_params),length(constraints));
Bic = struct(length(nb_params),length(constraints));

for c1 = 1:length(constraints)
    for c2 = 1:length(nb_params)

        % Choose number of parameters and constrint on covariance matrices
        nb_param = nb_params(c2);
        Model.cstr.Sigma = constraints{c1};

        % Generate training dataset
        [Xtrain, Ytrain] = GenerateScalableSignals(p(1:nb_param), int, nb_train_signals, 'qRandom');
        [Xtest,  Ytest]  = GenerateScalableSignals(p(1:nb_param), int, nb_test_signals, 'Random');
        
        % Add noise
        Xtrain  = AddNoise(Xtrain, Model.snrtrain);
        Xtest   = AddNoise(Xtest, 50);

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
                [theta,~,ll] = EstimateInverseFunction(Ytrain, Xtrain, K(k), Model.Lw, Model.maxiter, Model.cstr, 0);
                Ypredict	= EstimateParametersFromModel(Xtest, theta, 0);

                Rmse(k,:)   = EvaluateEstimation(Ytest, Ypredict);

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
                warning(['Iteration ' Model.cstr.Sigma '/' num2str(nb_param) '/' num2str(K(k)) ' failed'])
            end
        end

        Err{c2,c1}  = mean(Rmse,2);
        Bic{c2,c1}  = bic;
    end
end


%% Backup

if backup == 1
    clear tmp* Dico* Estim X* Y*
    save(['temp/' 'Kimpact'])
end


%% Disp

fig = figure;
hold on
for c1 = 1:length(constraints)
    subplot(1,length(constraints),c1)
    
    for c2 = 1:length(nb_params)
        plot(K, mean(Err{c2,c1},2), '.-', 'linewidth',1.5);
        %plot(K, mean(Dic{c2,c1},2), '.-', 'linewidth',1.5);
    end
end
ylabel('Mean RMSE')
%ylabel('BIC')
xlabel('K')
legend([repelem('P = ',length(nb_params),1) num2str(nb_params')])


%%

if backup == 1
    savefig(fig, ['figures/' 'Kimpact'])
end