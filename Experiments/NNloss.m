
%% Description
%
% Explaination
%
% Fabien Boux - 12/2020

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
nb_param = [3 4 5 6 7];

% Experiment settings
nb_test_signals = 1e4;
nb_train_signals = 1e4;
nb_params = [3 4 5 7];

% Regression settings
H = 6; 
Z = 100;
NeuralNet_.Exec = 'cpu';


%% Creating data

% Adding to path
addpath(genpath('functions'))
addpath(genpath('tools'))


%% Processing 

for n = 1:length(nb_params)
    
    % Generate signals
    [Xtrain, Ytrain] = GenerateScalableSignals(p(1:nb_param(n)), int, nb_train_signals, 'qRandom');
    [Xtest,  Ytest]  = GenerateScalableSignals(p(1:nb_param(n)), int, nb_test_signals, 'Random');
    
    %Neural network architecture
    layers = [sequenceInputLayer(size(Xtrain,2))];
    for h = 1:H
        layers = [layers 
            fullyConnectedLayer(Z)
            reluLayer; %tanhLayer 
        ];
    end
    layers = [layers
        fullyConnectedLayer(size(Ytrain,2))
        regressionLayer
        ];

    %Training settings
    options = trainingOptions('adam', ...
        'InitialLearnRate',1e-3, ...
        'SquaredGradientDecayFactor',0.9, ...
        'MaxEpochs',2000,...
        'MiniBatchSize',16, ...
        'Plots','none',...'training-progress', ...
        'ExecutionEnvironment', NeuralNet_.Exec,...
        'ValidationData',{Xtest',Ytest'},...
        'Verbose',0); %false

    %Train
    [~,info] = trainNetwork(Xtrain', Ytrain', layers, options);
    
    loss(n,:) = info.TrainingLoss;
    rmse(n,:) = info.TrainingRMSE;
end


%% Saving 

if backup == 1
    clear tmp* Xtrain Ytrain
    save(['temp/' mfilename])
end


%% Displaying

fig = figure;
% subplot(121)
hold on
plot(rmse', '.-')
ylabel('RMSE')
xlabel('Epoch')
% subplot(122)
% hold on
% plot(loss')


%% Exporting figures

if backup == 1
    savefig(fig, ['figures/' 'NNloss'])
end

