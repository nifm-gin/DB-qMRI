
%% Description
%
% We investigate the impact of adding noise on dictionary signals on the
% estimates of the DBL method for different dictionary sizes and SNR.
%
% Fabien Boux - 01/2020

Init
disp(['Running experiment ' mfilename '.m'])


%% Setting

% Execution settings
verbose = 1;
backup  = 1;

% Signal settings
int   	= [0.01 1];
p       = [.01 .01];
while min(abs(pdist(p',@(x,y) x-y))) < .05, p = 0.1 + 0.9*rand(1,10); end
nb_param = [3 5 7];
nb_signals = 4:5;
nb_signals = repmat(nb_signals,length(nb_param),1);
for i = 1:length(nb_param), nb_signals(i,:) = nb_signals(i,:).^nb_param(i); end

% Experiment settings
snr_on_dico = [inf 90 60 30];
snr_levels = [logspace(1, 2.035, 19) inf];
nb_test_signals = 10000;

% Regression settings
% db_method = 'DB-SL';
% Model_.K = 50;
db_method = 'DB-DL';
Model_.Exec = 'cpu';


%% Creating data

mNRMSE1 = nan(length(snr_levels),size(nb_signals,1),size(nb_signals,2));
                   mRMSE1 = mNRMSE1; mMAE1 = mNRMSE1;
mNRMSE2 = mNRMSE1; mRMSE2 = mNRMSE1; mMAE2 = mNRMSE1;
mNRMSE3 = mNRMSE1; mRMSE3 = mNRMSE1; mMAE3 = mNRMSE1;
mNRMSE4 = mNRMSE1; mRMSE4 = mNRMSE1; mMAE4 = mNRMSE1;
mNRMSE5 = mNRMSE1; mRMSE5 = mNRMSE1; mMAE5 = mNRMSE1;


%% Processing

for n = 1:size(nb_signals,1)
    
    for f = 1:size(nb_signals,2)
        
        if verbose >= 1, disp([num2str(nb_param(n)) ' - ' num2str(nb_signals(n,f))]); end
        
        % Generate training dataset
        [X, Y]  = GenerateScalableSignals(p(1:nb_param(n)), int, nb_signals(n,f), 'qRandom');
        
        %Learning
        Dico    = FormatDico(AddNoise(X, snr_on_dico(1)), Y);
        [~,Model1] = AnalyzeMRImages([], Dico,db_method, Model_);
        
        Dico    = FormatDico(AddNoise(X, snr_on_dico(2)), Y);
        [~,Model2] = AnalyzeMRImages([], Dico,db_method, Model_);
        
        Dico    = FormatDico(AddNoise(X, snr_on_dico(3)), Y);
        [~,Model3] = AnalyzeMRImages([], Dico,db_method, Model_);
        
        Dico    = FormatDico(AddNoise(X, snr_on_dico(4)), Y);
        [~,Model4] = AnalyzeMRImages([], Dico,db_method, Model_);
        
        Dico    = FormatDico(AddNoise(X, 10+90*rand(size(X,1),1)), Y);
        [~,Model5] = AnalyzeMRImages([], Dico,db_method, Model_);
        
        
        %Estimation
        parfor snr = 1:length(snr_levels)

            if verbose == 2, disp(['\t (n,f) = (' num2str(n) ',' num2str(f) ')\t Snr order: ' num2str(snr_levels(snr))]); end
            
            % Generate test data
            [Xtest, Ytest] = GenerateScalableSignals(p(1:nb_param(n)), int, nb_test_signals, 'Random');
            
            % Add noise on 
            [XtestN, tmp]       = AddNoise(Xtest, snr_levels(snr));
            real_snr(snr,n,f)   = mean(tmp);
            
            % Compute regressions
            % No noise
            Estim               = AnalyzeMRImages(XtestN, [], db_method,Model1, Ytest(:,1:size(Dico.Parameters.Par,2)));
            mNRMSE1(snr,n,f)    = mean(Estim.Regression.Errors.Nrmse);
            mRMSE1(snr,n,f)     = mean(Estim.Regression.Errors.Rmse);
            mMAE1(snr,n,f)      = mean(Estim.Regression.Errors.Mae);
            
            % 1st SNR
            Estim               = AnalyzeMRImages(XtestN, [], db_method,Model2, Ytest(:,1:size(Dico.Parameters.Par,2)));
            mNRMSE2(snr,n,f)    = mean(Estim.Regression.Errors.Nrmse);
            mRMSE2(snr,n,f)     = mean(Estim.Regression.Errors.Rmse);
            mMAE2(snr,n,f)      = mean(Estim.Regression.Errors.Mae);
            
            % 2nd SNR
            Estim               = AnalyzeMRImages(XtestN, [], db_method,Model3, Ytest(:,1:size(Dico.Parameters.Par,2)));
            mNRMSE3(snr,n,f)    = mean(Estim.Regression.Errors.Nrmse);
            mRMSE3(snr,n,f)     = mean(Estim.Regression.Errors.Rmse);
            mMAE3(snr,n,f)      = mean(Estim.Regression.Errors.Mae);
            
            % 3rd SNR
            Estim               = AnalyzeMRImages(XtestN, [], db_method,Model4, Ytest(:,1:size(Dico.Parameters.Par,2)));
            mNRMSE4(snr,n,f)    = mean(Estim.Regression.Errors.Nrmse);
            mRMSE4(snr,n,f)     = mean(Estim.Regression.Errors.Rmse);
            mMAE4(snr,n,f)      = mean(Estim.Regression.Errors.Mae);
            
            % range SNR
            Estim               = AnalyzeMRImages(XtestN, [], db_method,Model5, Ytest(:,1:size(Dico.Parameters.Par,2)));
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

for n = 2:size(nb_signals,1)
    for f = 1:size(nb_signals,2)
        
        ax(f) = subplot(size(nb_signals,1)-1,size(nb_signals,2), size(nb_signals,2) *(n-1-1) + f);
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
    savefig(fig, ['figures/' 'NoisyDicoSignals'])
end
