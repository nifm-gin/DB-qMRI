
%% Description
%
% We investigate the impact of dictionary size on the estimation accuracy,
% we generated dictionaries with different sizes and considered different
% number of parameters. We also characterize the impact of noise and we
% measure the computation time.
%
% Fabien Boux - 01/2020


%% Setting

% Execution settings
verbose = 1; %0, 1 or 2 for more details
backup  = 1;

% Signal settings
int     = [0 1];
p       = [.01 .01];
while min(abs(pdist(p',@(x,y) x-y))) < .05, p = 0.1 + 0.9*rand(1,10); end
nb_param = [3 4 5 6 7];
nb_signals = 3:6;
nb_signals = repmat(nb_signals,length(nb_param),1);
for i = 1:length(nb_param), nb_signals(i,:) = nb_signals(i,:).^nb_param(i); end

% Experiment settings
snr_levels = [logspace(1, 2.035, 40) inf];
nb_test_signals = 10000;

% Regression settings
Parameters = [];
Parameters.K = 50;
Parameters.cstr.Sigma  = 'd*';
Parameters.cstr.Gammat = ''; 
Parameters.cstr.Gammaw = '';
Parameters.Lw = 0;
snr_train = 60;

fast_limit = 500;


%% Creating data

% Adding to path
addpath(genpath('functions'))
addpath(genpath('tools'))

% Creating 


%% Processing 

mNRMSE_grid  = nan(length(snr_levels), size(nb_signals,1), size(nb_signals,2));
mRMSE_grid   = mNRMSE_grid; mMAE_grid = mNRMSE_grid; t_grid = mNRMSE_grid;
mNRMSE_gllim = nan(length(snr_levels), size(nb_signals,1), size(nb_signals,2));
mRMSE_gllim  = mNRMSE_gllim; mMAE_gllim = mNRMSE_gllim; t_gllim = mNRMSE_gllim;
Steps        = nan(size(nb_signals,1), size(nb_signals,2));

if verbose >= 1, disp('Nber of parameters - Dictionary size'); end

for n = 1:size(nb_signals,1)
    
    % Prepare outliers structure
    clear outliers
    for r = 1:nb_param(n), outliers{r} = int; end

    for f = 1:size(nb_signals,2)
        
        if verbose >= 1, disp([num2str(nb_param(n)) ' - ' num2str(nb_signals(n,f))]); end
        
        % Compute dico grid
        clear X Y
        if nb_param(n) ~= 1
            nb_step = nb_signals(n,f)^(1/nb_param(n));
        else
            nb_step = nb_signals(n,f);
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
        DicoR{1}.MRSignals = AddNoise(abs(X), snr_train); 
        DicoR{1}.Parameters.Par = Y;

        if size(DicoG{1}.MRSignals,1) ~= size(DicoR{1}.MRSignals,1)
            warning('Sizes are not equals')
        end
        
        if nb_signals(n,f) > fast_limit
            Parameters.K = 50;
        else
            Parameters.K = 10;
        end
                
        parfor snr = 1:length(snr_levels)

            if verbose == 2, disp(['\t (n,f) = (' num2str(n) ',' num2str(f) ')\t Snr order: ' num2str(snr_levels(snr))]); end
            
            % Generate test data
            Ytest 	= int(1) + (int(2) - int(1)) * rand(nb_test_signals,nb_param(n));
            Xtest = [];
            for sim = 1:size(Ytest,1)
                Xtest(sim,:) = toyMRsignal(Ytest(sim,:),p(1:nb_param(n)));
            end
            
            % Add noise
            [XtestN, tmp]       = AddNoise(Xtest, snr_levels(snr));
            real_snr(snr,n,f)   = mean(tmp); tmp = [];

            % Perform DBM
            tic;
            Estim   = AnalyzeMRImages(XtestN,DicoG,'DBM',[],Ytest(:,1:size(DicoG{1}.Parameters.Par,2)));
            
            t_grid(snr,n,f)         = toc;
            mNRMSE_grid(snr,n,f)    = nanmean(Estim.GridSearch.Errors.Nrmse);
            mRMSE_grid(snr,n,f)     = nanmean(Estim.GridSearch.Errors.Rmse);
            mMAE_grid(snr,n,f)      = nanmean(Estim.GridSearch.Errors.Mae);
            
            % Perform DBL
            tic;
            Estim   = AnalyzeMRImages(XtestN,DicoR,'DBL',Parameters,Ytest(:,1:size(DicoR{1}.Parameters.Par,2)),outliers);
            
            t_gllim(snr,n,f)        = toc;
            mNRMSE_gllim(snr,n,f)   = nanmean(Estim.Regression.Errors.Nrmse);
            mRMSE_gllim(snr,n,f)    = nanmean(Estim.Regression.Errors.Rmse);
            mMAE_gllim(snr,n,f)     = nanmean(Estim.Regression.Errors.Mae);
            
        end %snr
    end
end

t(:,1,:,:)      = t_grid;
t(:,2,:,:)      = t_gllim;
mNRMSE(:,1,:,:) = mNRMSE_grid;
mNRMSE(:,2,:,:) = mNRMSE_gllim;
mRMSE(:,1,:,:)  = mRMSE_grid;
mRMSE(:,2,:,:)  = mRMSE_gllim;
mMAE(:,1,:,:)   = mMAE_grid;
mMAE(:,2,:,:)   = mMAE_gllim;


%% Saving 

if backup == 1
    clear tmp* Dico* Estim
    save(['temp/' 'DictionarySize'])
end


%% Displaying

nb_param_wted_v = [5 7];

colors  = [          0    0.4470    0.7410
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
        plot([0 200], [mRMSE(end,1,nb_param_wted_,f) mRMSE(end,1,nb_param_wted_,f)], '--', 'linewidth',1.5, 'Color',colors(f,:))     
        pp(f) = plot(real_snr(:,nb_param_wted_,f), mRMSE(:,1,nb_param_wted_,f), '.-', 'MarkerSize', 18, 'Color', colors(f,:));          
        
        axis tight;
        
        ylim([0 .25]);
        xlim([11 105])
    end
    set(gca,'FontSize',16)
    switch count
        case 1
            title('DBM \newline (a)');
        case 3
            title('(c)');           
    end
    if count > 2, xlabel('SNR'); end
    ylabel([num2str(nb_param_wted) ' parameters\newline Average RMSE (s)'])
    
    ax(count+1) = subplot(2,2,count+1);
    hold on
    for f = 1:size(nb_signals,2)
        plot([0 200], [mRMSE(end,2,nb_param_wted_,f) mRMSE(end,2,nb_param_wted_,f)], '--', 'linewidth',1.5, 'Color',colors(f,:),'HandleVisibility','off')        
        pp(f) = plot(real_snr(:,nb_param_wted_,f), mRMSE(:,2,nb_param_wted_,f), '.-', 'MarkerSize', 18, 'Color',colors(f,:));
        
        axis tight;
        ylim([0 .25]);
        xlim([11 105])
    end
    
    if count > 2, xlabel('SNR'); end
    switch count
        case 1
            title('DBL \newline (b)');
        case 3
            title('(d)');           
    end
    count = 3;
    
    set(gca,'FontSize',16)
    lgd.FontSize = 16;
    
    legend([num2str(nb_signals(nb_param_wted_,:)') repelem('-signal dictionary',length(nb_signals(nb_param_wted_,:)),1)])
end

linkaxes(ax(1:2), 'xy')
linkaxes(ax(3:4), 'xy')
% set(ax,'XScale','log')


%% Supplementary figure

fig_supp = figure;

ccount = 1;
for n = 1:size(nb_signals,1)
    
    count = 1;
    for f = 1:size(nb_signals,2)
        
        if nb_signals(n,f) >= 100 % If we have at least 100 signals 
            
            ax(ccount) = subplot(size(nb_signals,1),size(nb_signals,2), size(nb_signals,2) *(n-1) + f);
            ccount= ccount +1;
            set(groot,'defaultAxesColorOrder',colors)
            
            hold on;
            plot([0 200], [mRMSE(end,1,n,f) mRMSE(end,1,n,f)], '--', 'linewidth',1.5, 'Color',colors(1,:),'HandleVisibility','off')        
            plot([0 200], [mRMSE(end,2,n,f) mRMSE(end,2,n,f)], '--', 'linewidth',1.5, 'Color',colors(2,:),'HandleVisibility','off')        
        
            plot(real_snr(:,n,f), mRMSE(:,1,n,f), '.-', 'MarkerSize', 12, 'Color', colors(1,:))
            plot(real_snr(:,n,f), mRMSE(:,2,n,f), '.-', 'MarkerSize', 12, 'Color', colors(2,:))
            
            axis tight;

            if f ==  size(nb_signals,2) && n == 1, legend('DBM','DBL'); end
            title([num2str(nb_signals(n,f)) ' signals']);
            if count == 1, ylabel([num2str(nb_param(n)) ' parameters \newline Average RMSE (s)']); end
            ylim([0 .21]);
            xlim([11 105])

            set(gca,'FontSize',12)
            lgd.FontSize = 16;
            
            count = count +1;
        end
        
        if n == size(nb_signals,1), xlabel('SNR'); end
    end
end
% set(ax,'XScale','log')


%% Exporting figures

if backup == 1
    savefig(fig, ['figures/' 'DictionarySize'])
    savefig(fig_supp, ['figures/' 'DictionarySize' '-supp'])
end

