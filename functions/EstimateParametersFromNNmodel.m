function [Ypredict] = EstimateParametersFromNNmodel(Xtest,net,gpu_opt)


if ~exist('gpu_opt','var')
    gpu_opt = 'auto';
end

%%

Ypredict = predict(net,Xtest','ExecutionEnvironment',gpu_opt)'; 

end

