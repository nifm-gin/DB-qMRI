
%% Description
%
% We investigate the impact of adding noise on dictionary signals on the
% estimates of the DBL method for different dictionary sizes and SNR.
%
% Fabien Boux - 01/2020


%% Setting

% Execution settings
verbose = 1;
backup  = 1;

% Signal settings
int   	= [0 1];
p       = [.01 .01];
while min(abs(pdist(p',@(x,y) x-y))) < .05, p = 0.1 + 0.9*rand(1,10); end
nb_param = [3 5 7];
nb_signals = 3:5;
nb_signals = repmat(nb_signals,length(nb_param),1);
for i = 1:length(nb_param), nb_signals(i,:) = nb_signals(i,:).^nb_param(i); end

% Experiment settings
snr_on_dico = [inf 90 60 30];
snr_levels = [logspace(1, 2.035, 19) inf];
nb_test_signals = 10000;

% Regression settings
Parameters.K = 50;
Parameters.cstr.Sigma  = 'd*';
Parameters.cstr.Gammat = ''; 
Parameters.cstr.Gammaw = '';
Parameters.Lw = 0;


%% Processing

mNRMSE1 = nan(length(snr_levels),size(nb_signals,1),size(nb_signals,2)); mRMSE1 = mNRMSE1; mMAE1 = mNRMSE1;
mNRMSE2 = nan(length(snr_levels),size(nb_signals,1),size(nb_signals,2)); mRMSE2 = mNRMSE1; mMAE2 = mNRMSE1;
mNRMSE3 = nan(length(snr_levels),size(nb_signals,1),size(nb_signals,2)); mRMSE3 = mNRMSE1; mMAE3 = mNRMSE1;
mNRMSE4 = nan(length(snr_levels),size(nb_signals,1),size(nb_signals,2)); mRMSE4 = mNRMSE1; mMAE4 = mNRMSE1;
mNRMSE5 = nan(length(snr_levels),size(nb_signals,1),size(nb_signals,2)); mRMSE5 = mNRMSE1; mMAE5 = mNRMSE1;

for n = 1:size(nb_signals,1)
    
    % Prepare outliers structure
    clear outliers
    for r = 1:nb_param(n), outliers{r} = int; end

    for f = 1:size(nb_signals,2)
        
        if verbose >= 1, disp([num2str(nb_param(n)) ' - ' num2str(nb_signals(n,f))]); end
        
        % Generate training dataset
        clear X Y
        Y       = int(1) + (int(2) - int(1)) * net(scramble(sobolset(nb_param(n)),'MatousekAffineOwen'),nb_signals(n,f));
        for sim = 1:size(Y,1)
            X(sim,:) = toyMRsignal(Y(sim,:),p(1:nb_param(n)));
        end
        DicoR{1}.MRSignals = abs(X); 
        DicoR{1}.Parameters.Par = Y;
        clear X Y        
        
        parfor snr = 1:length(snr_levels)

            if verbose == 2, disp(['\t (n,f) = (' num2str(n) ',' num2str(f) ')\t Snr order: ' num2str(snr_levels(snr))]); end
            
            % Generate test signals
            Ytest       = int(1) + (int(2) - int(1)) * rand(nb_test_signals,nb_param(n));

            Xtest = [];
            for sim = 1:size(Ytest,1)
                Xtest(sim,:) = toyMRsignal(Ytest(sim,:),p(1:nb_param(n)));
            end
            
            % Add noise on 
            [XtestN, tmp]       = AddNoise(Xtest, snr_levels(snr));
            real_snr(snr,n,f)   = mean(tmp);
            
            % Compute regressions
            % No noise
            Dico = DicoR;
            Dico{1}.MRSignals   = AddNoise(DicoR{1}.MRSignals,snr_on_dico(1));
            Estim               = AnalyzeMRImages(XtestN,Dico,'DBL',Parameters,Ytest(:,1:size(Dico{1}.Parameters.Par,2)),outliers);
            mNRMSE1(snr,n,f)  = mean(Estim.Regression.Errors.Nrmse);
            mRMSE1(snr,n,f)   = mean(Estim.Regression.Errors.Rmse);
            mMAE1(snr,n,f)    = mean(Estim.Regression.Errors.Mae);
            
            % 1st SNR
            Dico{1}.MRSignals   = AddNoise(DicoR{1}.MRSignals,snr_on_dico(2));
            Estim               = AnalyzeMRImages(XtestN,Dico,'DBL',Parameters,Ytest(:,1:size(Dico{1}.Parameters.Par,2)),outliers);
            mNRMSE2(snr,n,f)    = mean(Estim.Regression.Errors.Nrmse);
            mRMSE2(snr,n,f)     = mean(Estim.Regression.Errors.Rmse);
            mMAE2(snr,n,f)      = mean(Estim.Regression.Errors.Mae);
            
            % 2nd SNR
            Dico{1}.MRSignals   = AddNoise(DicoR{1}.MRSignals,snr_on_dico(3));
            Estim               = AnalyzeMRImages(XtestN,Dico,'DBL',Parameters,Ytest(:,1:size(Dico{1}.Parameters.Par,2)),outliers);
            mNRMSE3(snr,n,f)    = mean(Estim.Regression.Errors.Nrmse);
            mRMSE3(snr,n,f)     = mean(Estim.Regression.Errors.Rmse);
            mMAE3(snr,n,f)      = mean(Estim.Regression.Errors.Mae);
            
            % 3rd SNR
            Dico{1}.MRSignals   = AddNoise(DicoR{1}.MRSignals,snr_on_dico(4));
            Estim               = AnalyzeMRImages(XtestN,Dico,'DBL',Parameters,Ytest(:,1:size(Dico{1}.Parameters.Par,2)),outliers);
            mNRMSE4(snr,n,f)    = mean(Estim.Regression.Errors.Nrmse);
            mRMSE4(snr,n,f)     = mean(Estim.Regression.Errors.Rmse);
            mMAE4(snr,n,f)      = mean(Estim.Regression.Errors.Mae);
            
            % range SNR
            Dico{1}.MRSignals   = AddNoise(DicoR{1}.MRSignals, 10+90*rand(size(DicoR{1}.MRSignals,1),1));
            Estim               = AnalyzeMRImages(XtestN,Dico,'DBL',Parameters,Ytest(:,1:size(Dico{1}.Parameters.Par,2)),outliers);
            mNRMSE5(snr,n,f)    = mean(Estim.Regression.Errors.Nrmse);
            mRMSE5(snr,n,f)     = mean(Estim.Regression.Errors.Rmse);
            mMAE5(snr,n,f)      = mean(Estim.Regression.Errors.Mae);
            
        end %snr
    end
end


%% Saving 

if backup == 1
    clear Dico* Estim X* Y*
    save(['temp/' 'NoisyDicoSignals'])
end


%% Displaying

colors = [           0    0.4470    0.7410
                0.8500    0.3250    0.0980
                0.9290    0.6940    0.1250
                0.4940    0.1840    0.5560
                0.4660    0.6740    0.1880
                0.3010    0.7450    0.9330
                0.6350    0.0780    0.1840];
          
fig = figure;

for n = 1:size(nb_signals,1)
    for f = 1:size(nb_signals,2)
        
        ax(f) = subplot(size(nb_signals,1),size(nb_signals,2), size(nb_signals,2) *(n-1) + f);
        hold on
        plot([0 120], [mRMSE1(end,n,f) mRMSE1(end,n,f)], '--', 'linewidth',1.5, 'color', colors(1,:),'HandleVisibility','off')  
        plot(real_snr(:,n,f), mRMSE1(:,n,f), '.-','markersize',15, 'color', colors(1,:))
        
        plot([0 120], [mRMSE2(end,n,f) mRMSE2(end,n,f)], '--', 'linewidth',1.5, 'color', colors(2,:),'HandleVisibility','off')  
        plot(real_snr(:,n,f), mRMSE2(:,n,f),'.-','markersize',15, 'color', colors(2,:))
        
        plot([0 120], [mRMSE3(end,n,f) mRMSE3(end,n,f)], '--', 'linewidth',1.5, 'color', colors(3,:),'HandleVisibility','off')  
        plot(real_snr(:,n,f), mRMSE3(:,n,f), '.-','markersize',15, 'color', colors(3,:))
        
        plot([0 120], [mRMSE4(end,n,f) mRMSE4(end,n,f)], '--', 'linewidth',1.5, 'color', colors(4,:),'HandleVisibility','off')  
        plot(real_snr(:,n,f), mRMSE4(:,n,f),'.-','markersize',15, 'color', colors(4,:))
        
        plot([0 120], [mRMSE5(end,n,f) mRMSE5(end,n,f)], '--', 'linewidth',1.5, 'color', colors(5,:),'HandleVisibility','off')  
        plot(real_snr(:,n,f), mRMSE5(:,n,f), '.-','markersize',15, 'color', colors(5,:))
         
        if f == 1
            ylabel([num2str(nb_param(n)) ' parameters \newline Average RMSE (s)'])
        end
        if f == size(nb_signals,2)
           legend([repelem('SNR_{dico} = ',4,1) num2str(snr_on_dico(1:end)'); 'SNR_{dico} = [ ]'])
        end
        title([num2str(nb_signals(n,f)) ' signals'])
        
        if n == size(nb_signals,1), xlabel('SNR'); end
        xlim([10 115])
    end
    linkaxes(ax, 'y')
    
%     set(ax,'XScale','log')
%     set(ax,'fontsize',14)
end


%% Exporting figures

if backup == 1
    savefig(fig, ['figures/' 'NoisyDicoSignals-supp'])
end
