clear

%%
% Function definition
int	= [0 1];
p = [.01 .01];
while min(abs(pdist(p',@(x,y) x-y))) < .05
    p 	= 0.1 + 0.9*rand(1,10);
end % dishonest  or not ?
%p       = 0.12:0.12:1;

% Parameters definition
nb_param    = [3 5]; %7
nb_signals  = 3:5; %6
nb_signals  = repmat(nb_signals,length(nb_param),1);
for i = 1:length(nb_param), nb_signals(i,:) = nb_signals(i,:).^nb_param(i); end

snr_levels  = [logspace(1, 2.035, 9) inf]; %30

nb_test_signals = 10000;

fast_limit = 500;


%% Start computation

mNRMSE_grid  = nan(length(snr_levels), size(nb_signals,1), size(nb_signals,2));
mRMSE_grid   = mNRMSE_grid; mMAE_grid = mNRMSE_grid; t_grid = mNRMSE_grid;
mNRMSE_gllim = nan(length(snr_levels), size(nb_signals,1), size(nb_signals,2));
mRMSE_gllim  = mNRMSE_gllim; mMAE_gllim = mNRMSE_gllim; t_gllim = mNRMSE_gllim;
mNRMSE_nn    = nan(length(snr_levels), size(nb_signals,1), size(nb_signals,2));
mRMSE_nn     = mNRMSE_gllim; mMAE_nn = mNRMSE_gllim; t_nn = mNRMSE_gllim;
Steps        = nan(size(nb_signals,1), size(nb_signals,2));

%%

for n = 1:size(nb_signals,1)
    
    
    % Prepare outliers structure
    clear outliers
    for r = 1:nb_param(n), outliers{r} = int; end

    for f = 1:size(nb_signals,2)
        
        if n == size(nb_signals,1) && f == size(nb_signals,2), break; end
        
        disp([n f])
        
        % Compute dico grid
        clear X Y
        if nb_param(n) ~= 1
            nb_step    = nb_signals(n,f)^(1/nb_param(n));
        else
            nb_step    = nb_signals(n,f);
        end
        step    = (int(2)-int(1))/nb_step;
        v       = int(1)+step/2:step:int(2)-step/2;
        Y       = arrangement(v,nb_param(n));
        
        % Save the step size
        Steps(n,f) = step;
        
        for sim = 1:size(Y,1)
            X(sim,:) = toyMRsignal(Y(sim,:),p(1:nb_param(n)));
        end
        DicoG{1}.MRSignals = abs(X); 
        DicoG{1}.Parameters.Par = Y;
        clear X Y
        
        % Compute training dataset
        Y  	= int(1) + (int(2) - int(1)) * net(scramble(sobolset(nb_param(n)),'MatousekAffineOwen'),nb_signals(n,f));
        
        for sim = 1:size(Y,1)
            X(sim,:) = toyMRsignal(Y(sim,:),p(1:nb_param(n)));
        end
        DicoR{1}.MRSignals = abs(X); 
        DicoR{1}.Parameters.Par = Y;

        if size(DicoG{1}.MRSignals,1) ~= size(DicoR{1}.MRSignals,1)
            warning('Sizes are not equals')
        end
        
        
        parfor snr = 1:length(snr_levels)

            fprintf(['\t (n,f) = (' num2str(n) ',' num2str(f) ')\t Snr order: ' num2str(snr_levels(snr))  '\n'])
            
            % Generate test data
            Ytest       = int(1) + (int(2) - int(1)) * rand(nb_test_signals,nb_param(n));

            Xtest = [];
            for sim = 1:size(Ytest,1)
                Xtest(sim,:) = toyMRsignal(Ytest(sim,:),p(1:nb_param(n)));
            end
            
            % Add noise
            [XtestN, tmp]       = AddNoise(Xtest, snr_levels(snr));
            real_snr(snr,n,f)   = mean(tmp); tmp = [];

            % Compute regression and grid search
            tic;
            Estim           = AnalyzeMRImages(XtestN,DicoG,'ClassicMRF',[],Ytest(:,1:size(DicoG{1}.Parameters.Par,2)));
            mNRMSE_grid(snr,n,f) = nanmean(Estim.GridSearch.Errors.Nrmse);
            mRMSE_grid(snr,n,f) = nanmean(Estim.GridSearch.Errors.Rmse);
            mMAE_grid(snr,n,f) = nanmean(Estim.GridSearch.Errors.Mae);
            t_grid(snr) = toc;
            
            Parameters = [];
            if nb_signals(n,f) > fast_limit
                Parameters.K = 50;
                Parameters.Lw = 0;
                Parameters.cstr.Sigma  = 'd*';
            else
                Parameters.K = 10;
                Parameters.Lw = 0;
                Parameters.cstr.Sigma  = 'd*';
            end
            
            tic;
            Estim  	= AnalyzeMRImages(XtestN,DicoR,'RegressionMRF',Parameters,Ytest(:,1:size(DicoR{1}.Parameters.Par,2)),outliers);
            mNRMSE_gllim(snr,n,f) = nanmean(Estim.Regression.Errors.Nrmse);
            mRMSE_gllim(snr,n,f) = nanmean(Estim.Regression.Errors.Rmse);
            mMAE_gllim(snr,n,f) = nanmean(Estim.Regression.Errors.Mae);
            t_gllim(snr) = toc;
            
            tic;
            mn      = min(DicoR{1}.Parameters.Par );
            mx      = max(DicoR{1}.Parameters.Par  - mn);
            YtrainNN = (DicoR{1}.Parameters.Par  - mn) ./ mx;
            NeuralNet = EstimateNNmodel(DicoR{1}.MRSignals,YtrainNN,0);
            Y       = EstimateParametersFromNNmodel(XtestN,NeuralNet);
            Y       = (Y .* mx) + mn;
            [tmp1, tmp2, tmp3] = EvaluateEstimation(Ytest,Y);
            mRMSE_nn(snr,n,f)   = nanmean(tmp1);
            mNRMSE_nn(snr,n,f)  = nanmean(tmp2);
            mMAE_nn(snr,n,f)    = nanmean(tmp3);
            t_nn(snr) = toc;
            
        end %snr
    end
end


mNRMSE(:,1,:,:) = mNRMSE_grid;
mNRMSE(:,2,:,:) = mNRMSE_gllim;
mNRMSE(:,3,:,:) = mNRMSE_nn;
mRMSE(:,1,:,:)  = mRMSE_grid;
mRMSE(:,2,:,:)  = mRMSE_gllim;
mRMSE(:,3,:,:)  = mRMSE_nn;
mMAE(:,1,:,:)   = mMAE_grid;
mMAE(:,2,:,:)   = mMAE_gllim;
mMAE(:,3,:,:)   = mMAE_nn;


%% 

clear DicoG DicoR ax f1 funcfit l mdl pp
pause(5)
% 
% namefile = ['outputs/VsToy_' datestr(now,'yyyy-mmm-dd_HH:MM:SS') '.mat'];
% save(namefile)
