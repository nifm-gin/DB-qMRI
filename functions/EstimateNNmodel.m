function [net] = EstimateNNmodel(Xtrain,Ytrain,verb)


%% Create network layers

% K = floor(2/3 * size(Xtrain,2)) + size(Ytrain,2);
K = 300;

layers = [
    
    sequenceInputLayer(size(Xtrain,2))
    
    fullyConnectedLayer(K)
    tanhLayer%sigmoidLayer
    
    fullyConnectedLayer(K)
    tanhLayer
    
    fullyConnectedLayer(size(Ytrain,2))
    sigmoidLayer
    
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
    'InitialLearnRate',1e-4, ...
    'SquaredGradientDecayFactor',0.99, ...
    'MaxEpochs',2000, ...
    'MiniBatchSize',512, ...
    'Plots','none',...'training-progress', ...
    'Verbose',verb); %false

net = trainNetwork(Xtrain',Ytrain',layers,options);


end

