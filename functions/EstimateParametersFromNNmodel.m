function [Ypredict] = EstimateParametersFromNNmodel(Xtest,net)

Ypredict = predict(net,Xtest')'; 

end

