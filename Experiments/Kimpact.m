

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
Parameters = [];
Parameters.cstr.Sigma  = '';
Parameters.cstr.Gammat = ''; 
Parameters.cstr.Gammaw = '';
Parameters.Lw = 0;
snr_train  	= 60;


%% Creating data

% Adding to path
addpath(genpath('functions'))
addpath(genpath('tools'))

% for c = 1:length(constraints)
for c = 1:length(nb_params)
    
    nb_param = nb_params(c);
    Parameters.cstr.Sigma = constraints{1};

    
    Xtrain = []; Ytrain = [];
    Ytrain  	= int(1) + (int(2) - int(1)) * net(scramble(sobolset(nb_param),'MatousekAffineOwen'),nb_train_signals);
    parfor sim = 1:size(Ytrain,1)
        Xtrain(sim,:) = toyMRsignal(Ytrain(sim,:),p(1:nb_param));
    end
    DicoR = [];
    DicoR{1}.MRSignals = AddNoise(abs(Xtrain), snr_train); 
    DicoR{1}.Parameters.Par = Ytrain;

    Ytest       = int(1) + (int(2) - int(1)) * rand(nb_test_signals,nb_param);
    Xtest = [];
    parfor sim = 1:size(Ytest,1)
        Xtest(sim,:) = toyMRsignal(Ytest(sim,:),p(1:nb_param));
    end
    Xtest = AddNoise(Xtest, 100);
    
    Rmse    = zeros(length(K), size(Ytrain,2));
    bic     = zeros(1,length(K));
    [N,Lt]  = size(Ytrain);
    D       = size(Xtrain,2);
    
    
    for k = 1:length(K)

        disp(k)
        
        try
            [theta,~,ll]	= EstimateInverseFunction(Ytrain, Xtrain, K(k), Parameters.Lw, 200, Parameters.cstr, 0);
            Ypredict	= EstimateParametersFromModel(Xtest, theta, 0);

            Rmse(k,:) = EvaluateEstimation(Ytest, Ypredict);

            Lw = Parameters.Lw;
            switch Parameters.cstr.Sigma
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
        end
    end
    
    Err{c}  = mean(Rmse,2);
    Bic{c}  = bic;
end


%% Disp

fig = figure;
% subplot(121)
hold on
%for c = 1:length(constraints)
for c = 1:length(nb_params)
    plot(K, mean(Err{c},2), '.-', 'linewidth',1.5);
end
ylabel('Mean RMSE')
xlabel('K')

% subplot(122)
% hold on
% % for c = 1:length(constraints)
% for c = 1:length(nb_params)
%     plot(K, Bic{c}, '.-', 'linewidth',1.5)
% end
% title('BIC')
% % legend(constraints')
legend([repelem('P = ',length(nb_params),1) num2str(nb_params')])


%%

if backup == 1
    savefig(fig, ['figures/' 'Kimpact'])
end

save('temp_Kimpact')