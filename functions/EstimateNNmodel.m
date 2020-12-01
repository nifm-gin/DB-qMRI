function [net] = EstimateNNmodel(Xtrain,Ytrain,verb,gpu_opt,Z,H)


if ~exist('gpu_opt','var')
    gpu_opt = 'auto';
end

if ~exist('Z','var') 
    % Z = floor(2/3 * size(Xtrain,2)) + size(Ytrain,2);
    Z = 100;
end
if ~exist('H','var') 
    H = 6;
end


%% Create network layers

layers = [sequenceInputLayer(size(Xtrain,2))];
    
for h = 1:H
    
    layers = [layers 
        
        fullyConnectedLayer(Z)
        reluLayer; %tanhLayer 
    ];
end
    
layers = [layers
    
    fullyConnectedLayer(size(Ytrain,2))
    %reluLayer; %sigmoidLayer
    
    regressionLayer
    ];


%% Train network
%input: 25-point trajectory of magnitude images %V
%layer neural network containing two 300×300 hidden layers %V
%ADAM stochastic gradient descent algorithm %V
%loss function (cost) defined as the mean square error
%learning rate set to 0.001 %V
%sigmoid and tanh activation functions of the first and last hidden layers respectively %~
%output: 2 parameters T1 and T2 %V

options = trainingOptions('adam', ...
    'InitialLearnRate',1e-3, ...
    'SquaredGradientDecayFactor',0.9, ...
    'MaxEpochs',2000,...
    'MiniBatchSize',16, ...
    'Plots','none',...'training-progress', ...
    'ExecutionEnvironment', gpu_opt,...
    'Verbose',verb); %false

[net,info] = trainNetwork(Xtrain',Ytrain',layers,options);


end

