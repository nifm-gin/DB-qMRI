function [Estimation] = qDBM(Sequences, Dico, References)

if nargin < 2, error('Not enought input arguments'); end
if nargin < 3, References = []; end

% Format MRI data in 3D or 4D matrices
[Sequences,s1,s2,t,slices] = SequencesSizes(Sequences);
if ~isempty(References) && length(size(References))== 2
    References = reshape(References, s1, s2, size(References,2));
end

tic
for s = 1:slices
    
    %Estimation of parameters
    Estimation.GridSearch.Y(:,:,:,s) = ...
        reshape(EstimateParametersFromGrid(reshape(Sequences(:,:,:,s),s1*s2,t), Dico.MRSignals, Dico.Parameters.Par), s1,s2, []);

    %Errors computation if a reference image is provided
    if ~isempty(References)
        [Estimation.GridSearch.Errors.Rmse(s,:), Estimation.GridSearch.Errors.Nrmse(s,:), Estimation.GridSearch.Errors.Mae(s,:), Estimation.GridSearch.Errors.Nmae(s,:)] = ...
            EvaluateEstimation(reshape(References(:,:,:,s),s1*s2,size(References,3)), reshape(Estimation.GridSearch.Y(:,:,:,s),s1*s2,size(References,3)));
    end
end

Estimation.GridSearch.quantification_time = toc; 