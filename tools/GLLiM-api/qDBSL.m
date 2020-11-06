function [Estimation, Model] = qDBSL(Sequences, Dico, References, Model, SNRmap)

% Performed dictionary-based statistical learning (DB-SL) estimation using
% a dictionary, see:
% [Boux, Bayesian inverse regression for vascular magnetic resonance
% fingerprinting, 2020]
%
% Fabien Boux - 11/2020

% need to remove the outliers parameter

narginchk(2, 5);
if ~exist('References','var'),  References  = []; end
if ~exist('Model','var'),       Model       = []; end
if ~exist('SNRmap','var'),      SNRmap      = []; end

% Format MRI data in 3D or 4D matrices
[Sequences,s1,s2,t,slices] = SequencesSizes(Sequences);
if ~isempty(References) && length(size(References))== 2
    References = reshape(References, s1, s2, size(References,2));
end

% This parameters is only used to enable the parameter data normalisation
% it can possibly be set to zero but this is not recommended 
normalization = 1;

% If no model has already be computed (ie theta), learn it
tic
if ~any(strcmp(fieldnames(Model),'theta'))

    %Normalize training data (zero mean and unit variance)
    if normalization == 1
        Model.factors.Ymean	= nanmean(Dico.Parameters.Par,1);
        Model.factors.Ystd = nanstd(Dico.Parameters.Par);
        Dico.Parameters.Par = (Dico.Parameters.Par - Model.factors.Ymean) ./ Model.factors.Ystd;
        Model.factors.normalization = 1;
        
    else
        Model.factors.Ymean	= 0;
        Model.factors.Ystd = 1;
        Model.factors.normalization = 0;
    end

    Xtrain = Dico.MRSignals;
    [~,Model] = EstimateParametersFromRegression(Xtrain, Xtrain, Dico.Parameters.Par, Dico.Parameters.Par, Model);

else
    Dico.MRSignals       = [];
    Dico.Parameters.Par  = [];
end
Estimation.Regression.learning_time = toc; 

% Use the model to quantify estimates
tic
for s = 1:slices
            
    %Estimation of parameters
    if ~isempty(SNRmap)
        warning('off')

        %Recompute Sigma_k of the GLLiM model (if SNRMap given)
        Yestim = nan(s1*s2,size(Model.theta.A,2));
        Cov = nan(size(Model.theta.A,2),size(Model.theta.A,2),s1*s2);
        Pik = nan(s1*s2,Model.K);
        Mahaldist = nan(s1*s2,1);

        Param_updated = Model;
        SNRmap(SNRmap < 2) = 2;
        var_noise = (max(Sequences(:,:,:,s),[],3) ./ SNRmap(:,:,s)).^2;

        Nb_model = 20;
        var_bounds = sort(1 ./ [2:50/(0.75*Nb_model):62 62:50/(0.25*Nb_model):122 inf].^2);
        for m = 1:length(var_bounds)-1
            if any(reshape( (var_noise >= var_bounds(m)) & (var_noise < var_bounds(m+1)) ,1,[]))
                var_noise_ = var_bounds(m) + (var_bounds(m+1) - var_bounds(m))/2;
                Param_updated.theta = updateSigma(Model.theta,var_noise_);

                Sequences_in = Sequences(:,:,:,s) .* ( (var_noise >= var_bounds(m)) & (var_noise < var_bounds(m+1)) ) ;
                Sequences_in(Sequences_in == 0) = nan;
                [Yestim_in,~,Cov_in,~,~,Mahaldist_in] = ...
                    EstimateParametersFromRegression(reshape(Sequences_in(:,:,:),s1*s2,t), Dico.MRSignals, Dico.Parameters.Par, [], Param_updated);

                Yestim = nansum(cat(3,Yestim,Yestim_in),3);
                Cov = nansum(cat(4,Cov,Cov_in),4);
                Mahaldist = nansum(cat(3,Mahaldist,Mahaldist_in),3);
                %Pik = nansum(cat(3,Pik,Pik_in),3);
                Pik = [];
            end
        end
        warning('on')
    else
        [Yestim,~,Cov,~,Pik,Mahaldist] = ...
            EstimateParametersFromRegression(reshape(Sequences(:,:,:,s),s1*s2,t), Dico.MRSignals, Dico.Parameters.Par, [], Model);
    end

    %Rescale
    if any(strcmp(fieldnames(Model),'factors')) && normalization == 1
        Yestim(:,1:end-Model.Lw) = (Yestim(:,1:end-Model.Lw) .* Model.factors.Ystd) + Model.factors.Ymean;
    end
    Estimation.Regression.Y(:,:,:,s)    = reshape(Yestim(:,1:end-Model.Lw), s1,s2,[]);

    Cov  	= reshape(Cov,size(Cov,1),size(Cov,2),s1,s2);
    Cov     = Cov(1:end-Model.Lw,1:end-Model.Lw,:,:); 
    for ss = 1:s1
        for sss = 1:s2
            if any(strcmp(fieldnames(Model),'factors')) && normalization == 1
                Estimation.Regression.Cov(ss,sss,:,s)   = diag(Cov(:,:,ss,sss))';
                Estimation.Regression.Cov(ss,sss,:,s) = squeeze(Estimation.Regression.Cov(ss,sss,:,s))' .* Model.factors.Ystd.^2;
            else
                Estimation.Regression.Cov(ss,sss,:,s)   = diag(Cov(:,:,ss,sss))';
            end
        end
    end

    Estimation.Regression.Pik = Pik;
    Estimation.Regression.Mahaldist = Mahaldist;

    %Remove outliers
%     if ~isempty(Outliers)
%         if size(Estimation.Regression.Y,3) == length(Outliers)
%             for o = 1:length(Outliers)
%                 tmp = Estimation.Regression.Y(:,:,o,s);
%                 tmp(tmp < Outliers{o}(1)) = nan;
%                 tmp(tmp > Outliers{o}(2)) = nan;                        
%                 Estimation.Regression.Y(:,:,o,s) = tmp;
%             end
%         end
%     end

    %Errors computation if a reference image is provided
    if ~isempty(References)
        [Estimation.Regression.Errors.Rmse(s,:), Estimation.Regression.Errors.Nrmse(s,:), Estimation.Regression.Errors.Mae(s,:), Estimation.Regression.Errors.Nmae(s,:)] = ...
            EvaluateEstimation(reshape(References(:,:,:,s),s1*s2,size(References,3)), reshape(Estimation.Regression.Y(:,:,:,s),s1*s2,size(References,3)));
    end
end
Estimation.Regression.quantification_time = toc;