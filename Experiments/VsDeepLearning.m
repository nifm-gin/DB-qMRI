
%% Description
%
% Blabla
%
% Fabien Boux - 07/2020

Init
disp(['Running experiment ' mfilename '.m'])


%% Setting

% Execution settings
verbose = 1;
backup  = 1;

% Function definition
int     = [.01 1];
p       = [.01 .01];
while min(abs(pdist(p',@(x,y) x-y))) < .05, p = 0.1 + 0.9*rand(1,10); end

% Parameters definition
nb_param = [3 4 5 6 7];
nb_signals = 3:7;
nb_signals = repmat(nb_signals,length(nb_param),1);
for i = 1:length(nb_param), nb_signals(i,:) = nb_signals(i,:).^nb_param(i); end

% Experiment settings
snr_levels = [logspace(1, 2.035, 29) inf]; %39
nb_test_signals = 10000;

fast_limit = 500;
snr_train = 60;


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
        
        disp([n f])
        
        % Compute dico grid
        clear X Y
        nb_step    = nb_signals(n,f)^(1/nb_param(n));
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
        DicoR{1}.MRSignals = AddNoise(abs(X), snr_train); 
        DicoR{1}.Parameters.Par = Y;

        if size(DicoG{1}.MRSignals,1) ~= size(DicoR{1}.MRSignals,1)
            warning('Sizes are not equals')
        end
        
        if nb_signals(n,f) < 100000
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
                Estim           = AnalyzeMRImages(XtestN,DicoG,'DBM',[],Ytest(:,1:size(DicoG{1}.Parameters.Par,2)));
                mNRMSE_grid(snr,n,f) = nanmean(Estim.GridSearch.Errors.Nrmse);
                mRMSE_grid(snr,n,f) = nanmean(Estim.GridSearch.Errors.Rmse);
                mMAE_grid(snr,n,f) = nanmean(Estim.GridSearch.Errors.Mae);
                t_grid(snr) = toc;

                Parameters = [];
                if nb_signals(n,f) > fast_limit
                    Parameters.K = 50;
                else
                    Parameters.K = 10;
                end
                Parameters.cstr.Sigma  = 'd*';
                Parameters.cstr.Gammat = ''; 
                Parameters.cstr.Gammaw = '';
                Parameters.Lw = 0;

                tic;
                Estim  	= AnalyzeMRImages(XtestN,DicoR,'DBL',Parameters,Ytest(:,1:size(DicoR{1}.Parameters.Par,2)),outliers);
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
        elseif  nb_signals(n,f) < 500000
            for snr = 1:length(snr_levels)
                
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
                Estim           = AnalyzeMRImages(XtestN,DicoG,'DBM',[],Ytest(:,1:size(DicoG{1}.Parameters.Par,2)));
                mNRMSE_grid(snr,n,f) = nanmean(Estim.GridSearch.Errors.Nrmse);
                mRMSE_grid(snr,n,f) = nanmean(Estim.GridSearch.Errors.Rmse);
                mMAE_grid(snr,n,f) = nanmean(Estim.GridSearch.Errors.Mae);
                t_grid(snr) = toc;

                Parameters = [];
                if nb_signals(n,f) > fast_limit
                    Parameters.K = 50;
                else
                    Parameters.K = 10;
                end
                Parameters.cstr.Sigma  = 'd*';
                Parameters.cstr.Gammat = ''; 
                Parameters.cstr.Gammaw = '';
                Parameters.Lw = 0;

                tic;
                Estim  	= AnalyzeMRImages(XtestN,DicoR,'DBL',Parameters,Ytest(:,1:size(DicoR{1}.Parameters.Par,2)),outliers);
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

if backup == 1
    clear tmp* X* Y* Dic ax f1 funcfit l mdl pp
    save(['temp/' 'VsDeepLearning'])
end


%%

nb_param_wted_v = [3 5];

colors = [           0    0.4470    0.7410
                0.8500    0.3250    0.0980
                0.9290    0.6940    0.1250
                0.4940    0.1840    0.5560
                0.4660    0.6740    0.1880
                0.3010    0.7450    0.9330
                0.6350    0.0780    0.1840];
colors = [colors; colors];
    
funcfit = @(a,x) a(1) + a(2) * exp(-a(3)* x);


fig = figure;
count = 1;
for nb_param_wted = nb_param_wted_v
    
    nb_param_wted_ = find(nb_param == nb_param_wted);
    
    ax(count) = subplot(2,2,count);
    set(groot,'defaultAxesColorOrder',colors)
    hold on
    
    for f = 1:size(nb_signals,2)    

        pp(f) = plot(real_snr(:,nb_param_wted_,f), mRMSE(:,3,nb_param_wted_,f), '.', 'MarkerSize', 25, 'Color', colors(f,:));        
        plot([0 200], [mRMSE(end,3,nb_param_wted_,f) mRMSE(end,3,nb_param_wted_,f)], '-.', 'linewidth',1.5, 'Color',colors(f,:))        
        
        axis tight;
        
        ylim([.001 .25]);
        xlim([11 105])

%         conv_val(1,nb_param_wted_,f) = c(1,:);
    end
    set(gca,'FontSize',16)
    switch count
        case 1
            title('(a)');
        case 3
            title('(c)');           
    end
    if count > 2, xlabel('SNR'); end
    ylabel('Mean RMSE')
    
    ax(count+1) = subplot(2,2,count+1);
    hold on
    for f = 1:size(nb_signals,2)   
        pp(f) = plot(real_snr(:,nb_param_wted_,f), mRMSE(:,2,nb_param_wted_,f), '.', 'MarkerSize', 25, 'Color',colors(f,:));
        plot([0 200], [mRMSE(end,2,nb_param_wted_,f) mRMSE(end,2,nb_param_wted_,f)], '-.', 'linewidth',1.5, 'Color',colors(f,:))        
        
        legds{f}   = [num2str(nb_signals(nb_param_wted_,f)) '-signal dictionary'];

        axis tight;
        
        ylim([.001 .25]);
        xlim([11 105]) %???

        %conv_val(2,nb_param_wted_,f) = c(1,:);
    end
    
    if count > 2, xlabel('SNR'); end
    
    switch count
        case 1
            title('(b)');
        case 3
            title('(d)');           
    end
    
    count = 3;
    
    set(gca,'FontSize',16)
    %lgd = legend(pp, legds);
    lgd.FontSize = 16;
end

linkaxes(ax(1:2), 'xy')
linkaxes(ax(3:4), 'xy')


%%
if backup == 1
    savefig(fig, 'figures/VsDeepLearning.fig')
end


%%

fig_supp = figure;
for n = 1:size(nb_signals,1)
    
    count = 1;
    
    for f = 1:size(nb_signals,2)
        
        if nb_signals(n,f) >= 100 % If we have at least 100 signals 
            
            ax = subplot(size(nb_signals,1),size(nb_signals,2), size(nb_signals,2) *(n-1) + f);
            set(groot,'defaultAxesColorOrder',colors)
%             plot(real_snr(:,n,f), mRMSE(:,:,n,f), '.', 'MarkerSize', 15)
            
            %hold on;
            for i = 1:size(mRMSE,2)
                try
                    mdl = fitnlm(real_snr(:,n,f), mRMSE(:,i,n,f), funcfit, [0.3 0.4 0.1]);
                    c(:,i) = mdl.Coefficients.Estimate;
                    %plot(real_snr(:,n,f), funcfit(c,real_snr(:,n,f)))
                    aerror(:,i,n,f) = repelem(std(( mRMSE(:,i,n,f) - funcfit(c(:,i),real_snr(:,n,f)))),size(real_snr(:,n,f),1),1);
                end
            end
            %[l,p] = boundedline(real_snr(:,n,f)',funcfit(c(:,1),real_snr(:,n,f))', aerror(:,1,n,f)', '--', real_snr(:,n,f)',funcfit(c(:,2),real_snr(:,n,f))', aerror(:,2,n,f)', '--', 'cmap', colors);
            hold on;
            plot(real_snr(:,n,f), mRMSE(:,1,n,f), '.-', 'MarkerSize', 12)
            plot(real_snr(:,n,f), mRMSE(:,2,n,f), '.-', 'MarkerSize', 12)
            plot(real_snr(:,n,f), mRMSE(:,3,n,f), '.-', 'MarkerSize', 12)
            
            plot([0 120], [mRMSE(end,1,n,f) mRMSE(end,1,n,f)], '--', 'linewidth',1.5, 'Color',colors(1,:))        
            plot([0 120], [mRMSE(end,2,n,f) mRMSE(end,2,n,f)], '--', 'linewidth',1.5, 'Color',colors(2,:))        
            plot([0 120], [mRMSE(end,3,n,f) mRMSE(end,3,n,f)], '--', 'linewidth',1.5, 'Color',colors(3,:))        
            axis tight;

            if f ==  size(nb_signals,2) && n == 1, legend('Grid Search','Regression','NN'); end
    %         xlabel('SNR'); ylabel('Mean NRMSE')
            title([num2str(nb_signals(n,f)) ' signals']);
%             if count == 1, ylabel([num2str(nb_param(n)) ' parameters'], 'fontweight','bold'); end
%             ylim([0 .15]);
            %xlim([real_snr(1) real_snr(end)])

            conv_val(:,n,f) = c(1,:);

            set(gca,'FontSize',12)
            lgd.FontSize = 16;
            
            count = count +1;
        end
        
        if count >= 3, ylabel('Mean RMSE'); end
        if n == size(nb_signals,1), xlabel('SNR'); end
    end
    
end


%%
if backup == 1
    savefig(fig_supp, 'figures/VsDeepLearning-supp.fig')
end
