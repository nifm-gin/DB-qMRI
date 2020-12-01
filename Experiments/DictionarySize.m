
%% Description
%
% We investigate the impact of dictionary size on the estimation accuracy,
% we generated dictionaries with different sizes and considered different
% number of parameters. We also characterize the impact of noise and we
% measure the computation time.
%
% Fabien Boux - 01/2020

Init
disp(['Running experiment ' mfilename '.m'])


%% Setting

% Execution settings
verbose = 2; %0, 1 or 2 for more details
backup  = 1;

% Signal settings
int     = [.01 1];
p       = [.01 .01];
while min(abs(pdist(p',@(x,y) x-y))) < .05, p = 0.1 + 0.9*rand(1,10); end
nb_param = [3 4 5 6 7];
nb_signals = 3:6;
nb_signals = repmat(nb_signals,length(nb_param),1);
for i = 1:length(nb_param), nb_signals(i,:) = nb_signals(i,:).^nb_param(i); end

% Experiment settings
snr_levels = [logspace(1, 2.035, 39) inf];
nb_test_signals = 10000;

% Regression settings
K       = 50; 
fast_limit = 500;
dbdl_computation = 1;
NeuralNet_.Exec = 'cpu';


%% Init data

mNRMSE_grid = nan(length(snr_levels), size(nb_signals,1), size(nb_signals,2));
mRMSE_grid  = mNRMSE_grid; mMAE_grid = mNRMSE_grid; 
mNRMSE_gllim = nan(length(snr_levels), size(nb_signals,1), size(nb_signals,2));
mRMSE_gllim = mNRMSE_gllim; mMAE_gllim = mNRMSE_gllim;
mNRMSE_nn   = nan(length(snr_levels), size(nb_signals,1), size(nb_signals,2));
mRMSE_nn  	= mNRMSE_nn; mMAE_nn = mNRMSE_nn; 
t_grid   	= mNRMSE_grid;
t_gllim  	= mNRMSE_gllim; t_gllim_learn = mNRMSE_gllim;
t_nn        = mNRMSE_gllim; t_nn_learn = mNRMSE_gllim;
mem_size_grid = t_grid;
mem_size_gllim = mem_size_grid; 
mem_size_nn = mem_size_grid;


%% Processing 

if verbose >= 1, disp('Nber of parameters - Dictionary size'); end

for n = 1:size(nb_signals,1)
    
    % Prepare outliers structure
    clear outliers
    for r = 1:nb_param(n), outliers{r} = int; end

    for f = 1:size(nb_signals,2)
        
        if verbose >= 1, disp([num2str(nb_param(n)) ' - ' num2str(nb_signals(n,f))]); end
        
        % Generate dico grid
        [X, Y] = GenerateScalableSignals(p(1:nb_param(n)), int, nb_signals(n,f), 'Grid');
        Dico_DBM = FormatDico(X, Y);
        
        % Generate training dataset
        [X, Y] = GenerateScalableSignals(p(1:nb_param(n)), int, nb_signals(n,f), 'qRandom');
        Dico_DBL = FormatDico(X, Y);
        
        if size(Dico_DBM.MRSignals,1) ~= size(Dico_DBL.MRSignals,1)
            warning('Sizes are not equals')
        end
        
        if nb_signals(n,f) > fast_limit
            Model_.K = K;
        elseif (nb_signals(n,f) > 100) && (nb_signals(n,f) <= fast_limit)
            Model_.K = 30;
        else
            Model_.K = 10;
        end
        
        %learn models
        [~,Model] = AnalyzeMRImages([],Dico_DBL,'DB-SL',Model_);
        
        if dbdl_computation == 1
            [~,NeuralNet] = AnalyzeMRImages([],Dico_DBL,'DB-DL',NeuralNet_);
        end
        
        % Need to remove parfor loop for memory requirement computation
        for snr = 1:length(snr_levels)

            if verbose == 2, disp(['(n,f) = (' num2str(n) ',' num2str(f) ') - Snr order: ' num2str(snr_levels(snr))]); end
            
            % Generate test data
            [Xtest, Ytest] = GenerateScalableSignals(p(1:nb_param(n)), int, nb_test_signals, 'Random');
            
            % Add noise
            [XtestN, tmp]       = AddNoise(Xtest, snr_levels(snr));
            real_snr(snr,n,f)   = mean(tmp); tmp = [];

            % Perform DBM
            Estim   = AnalyzeMRImages(XtestN,Dico_DBM,'DBM',[],Ytest(:,1:size(Dico_DBM.Parameters.Par,2)));
            
            %time
            t_grid(snr,n,f)         = Estim.GridSearch.quantification_time;
                   
            %estimation accuracy
            mNRMSE_grid(snr,n,f)    = nanmean(Estim.GridSearch.Errors.Nrmse);
            mRMSE_grid(snr,n,f)     = nanmean(Estim.GridSearch.Errors.Rmse);
            mMAE_grid(snr,n,f)      = nanmean(Estim.GridSearch.Errors.Mae);
            
            % Perform DBL
            Estim   = AnalyzeMRImages(XtestN,[],'DB-SL',Model,Ytest(:,1:size(Dico_DBL.Parameters.Par,2)),outliers);

            %times
            t_gllim(snr,n,f)        = Estim.Regression.quantification_time;
            t_gllim_learn(snr,n,f)  = Estim.Regression.learning_time;
            
            %estimation accuracy
            mNRMSE_gllim(snr,n,f)   = nanmean(Estim.Regression.Errors.Nrmse);
            mRMSE_gllim(snr,n,f)    = nanmean(Estim.Regression.Errors.Rmse);
            mMAE_gllim(snr,n,f)     = nanmean(Estim.Regression.Errors.Mae);
            
            
            % Perform DB-DL
            if dbdl_computation == 1
                Estim   = AnalyzeMRImages(XtestN,[],'DB-DL',NeuralNet,Ytest(:,1:size(Dico_DBL.Parameters.Par,2)),outliers);

                %times
                t_nn(snr,n,f)        = Estim.Regression.quantification_time;
                t_nn_learn(snr,n,f)  = Estim.Regression.learning_time;
                
                %estimation accuracy
                mRMSE_nn(snr,n,f) 	= nanmean(Estim.Regression.Errors.Rmse);
                mNRMSE_nn(snr,n,f)  = nanmean(Estim.Regression.Errors.Nrmse);
                mMAE_nn(snr,n,f)    = nanmean(Estim.Regression.Errors.Mae);
            end
            
            %memory requirement
            me      = whos('Dico_DBM');
            mem_size_grid(snr,n,f)  = me.bytes;

            me      = whos('Model');
            mem_size_gllim(snr,n,f) = me.bytes;

            me      = whos('NeuralNet');
            mem_size_nn(snr,n,f)  = me.bytes;
            
        end %snr
    end
end

%times (sec)
t(:,1,:,:)      = t_grid;
t(:,2,:,:)      = t_gllim;
t(:,3,:,:)      = t_gllim_learn;

%memory (bits)
mem(:,1,:,:)  	= mem_size_grid;
mem(:,2,:,:)  	= mem_size_gllim;
mem(:,3,:,:)  	= mem_size_nn;

%errors
mNRMSE(:,1,:,:) = mNRMSE_grid;
mNRMSE(:,2,:,:) = mNRMSE_gllim;
mNRMSE(:,3,:,:) = mNRMSE_nn;
mRMSE(:,1,:,:)  = mRMSE_grid;
mRMSE(:,2,:,:)  = mRMSE_gllim;
mRMSE(:,3,:,:)  = mRMSE_nn;
mMAE(:,1,:,:)   = mMAE_grid;
mMAE(:,2,:,:)   = mMAE_gllim;
mMAE(:,3,:,:)   = mMAE_nn;


%% Saving 

if backup == 1
    clear tmp* Dico* Estim X* Y*
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
    
    ax(count) = subplot(2,3,count);
    set(groot,'defaultAxesColorOrder',colors)
    hold on
    
    for f = 1:size(nb_signals,2) 
        plot([1e-6 200], [mRMSE(end,1,nb_param_wted_,f) mRMSE(end,1,nb_param_wted_,f)], '--', 'linewidth',1.5, 'Color',colors(f,:))     
        pp(f) = plot(real_snr(:,nb_param_wted_,f), mRMSE(:,1,nb_param_wted_,f), '.-', 'MarkerSize', 18, 'Color', colors(f,:));          
        
        axis tight;
        
        ylim([0 .25]);
        xlim([11 105])
    end
    set(gca,'FontSize',16)
    switch count
        case 1
            title('DBM \newline (a)');
        case 4
            title('(d)');           
    end
    if count > 2, xlabel('SNR'); end
    ylabel([num2str(nb_param_wted) ' parameters\newline Average RMSE (s)'])
    
    
    ax(count+1) = subplot(2,3,count+1);
    hold on
    for f = 1:size(nb_signals,2)
        plot([1e-6 200], [mRMSE(end,2,nb_param_wted_,f) mRMSE(end,2,nb_param_wted_,f)], '--', 'linewidth',1.5, 'Color',colors(f,:),'HandleVisibility','off')        
        pp(f) = plot(real_snr(:,nb_param_wted_,f), mRMSE(:,2,nb_param_wted_,f), '.-', 'MarkerSize', 18, 'Color',colors(f,:));
        
        axis tight;
        ylim([0 .25]);
        xlim([11 105])
    end
    
    if count > 2, xlabel('SNR'); end
    switch count
        case 1
            title('DB-SL \newline (b)');
        case 4
            title('(e)');
    end
    set(gca,'FontSize',16)
    lgd.FontSize = 16;
    
    
    ax(count+2) = subplot(2,3,count+2);
    hold on
    for f = 1:size(nb_signals,2)
        plot([1e-6 200], [mRMSE(end,3,nb_param_wted_,f) mRMSE(end,3,nb_param_wted_,f)], '--', 'linewidth',1.5, 'Color',colors(f,:),'HandleVisibility','off')        
        pp(f) = plot(real_snr(:,nb_param_wted_,f), mRMSE(:,3,nb_param_wted_,f), '.-', 'MarkerSize', 18, 'Color',colors(f,:));
        
        axis tight;
        ylim([0 .25]);
        xlim([11 105])
    end
    
    if count > 2, xlabel('SNR'); end
    switch count
        case 1
            title('DB-DL \newline (c)');
        case 4
            title('(f)');
    end
    count = 4;
    
    set(gca,'FontSize',16)
    lgd.FontSize = 16;
    
    legend([num2str(nb_signals(nb_param_wted_,:)') repelem('-signal dictionary',length(nb_signals(nb_param_wted_,:)),1)])
end

linkaxes(ax(1:3), 'xy')
linkaxes(ax(4:5), 'xy')
%set(ax,'XScale','log')


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
            plot([1e-6 200], [mRMSE(end,1,n,f) mRMSE(end,1,n,f)], '--', 'linewidth',1.5, 'Color',colors(1,:),'HandleVisibility','off')        
            plot([1e-6 200], [mRMSE(end,2,n,f) mRMSE(end,2,n,f)], '--', 'linewidth',1.5, 'Color',colors(2,:),'HandleVisibility','off')        
            if dbdl_computation == 1
                plot([1e-6 200], [mRMSE(end,3,n,f) mRMSE(end,3,n,f)], '--', 'linewidth',1.5, 'Color',colors(3,:),'HandleVisibility','off')        
            end
            
            plot(real_snr(:,n,f), mRMSE(:,1,n,f), '.-', 'MarkerSize', 12, 'Color', colors(1,:))
            plot(real_snr(:,n,f), mRMSE(:,2,n,f), '.-', 'MarkerSize', 12, 'Color', colors(2,:))
            if dbdl_computation == 1
                plot(real_snr(:,n,f), mRMSE(:,3,n,f), '.-', 'MarkerSize', 12, 'Color', colors(3,:))
            end
            
            axis tight;

            if f ==  size(nb_signals,2) && n == 1
                if dbdl_computation == 1
                    legend('DBM','DB-SL','DB-DL');
                else
                    legend('DBM','DB-SL');
                end
            end
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
%set(ax,'XScale','log')


%% Exporting figures

if backup == 1
    savefig(fig, ['figures/' 'DictionarySize'])
    savefig(fig_supp, ['figures/' 'DictionarySize' '-supp'])
end


%% Main results

disp(' ')
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
 
% disp(['For 5 parameters the DBL RMSE is ' num2str(...
%     100* mean(reshape((mRMSE(:,1,3,:) - mRMSE(:,2,3,:))  ./ mRMSE(:,1,3,:), 1,[])), 3)...
%     ' +/- '  num2str(...
%     100* std(reshape((mRMSE(:,1,3,:) - mRMSE(:,2,3,:)) ./ mRMSE(:,1,3,:), 1,[])), 3)...
%     ' lower than the DBM RMSE'])
% disp(['For 7 parameters the DBL RMSE is ' num2str(...
%     100* mean(reshape((mRMSE(:,1,5,:) - mRMSE(:,2,5,:))  ./ mRMSE(:,1,5,:), 1,[])), 3)...
%     ' +/- '  num2str(...
%     100* std(reshape((mRMSE(:,1,5,:) - mRMSE(:,2,5,:)) ./ mRMSE(:,1,5   ,:), 1,[])), 3)...
%     ' lower than the DBM RMSE'])

w = diff(snr_levels);
w = w(1:end-1);
w = w' ./ sum(w);


disp(' ')
disp('Average RMSE obtained with the DBL method are lower than the one obtained')
disp('with the DBM method:')
disp([num2str(...
100* wmean(reshape((squeeze(mNRMSE(1:end-2,1,[3],:)) - squeeze(mNRMSE(1:end-2,2,[3],:)))...
    ./ squeeze(mNRMSE(1:end-2,1,[3],:)),1,[]), repmat(w,4,1)')...
, 3) ' +/- ' num2str(...
100* var(reshape((squeeze(mNRMSE(1:end-2,1,[3],:)) - squeeze(mNRMSE(1:end-2,2,[3],:)))...
    ./ squeeze(mNRMSE(1:end-2,1,[3],:)),1,[]), repmat(w,4,1)').^.5 ...
, 3) ' % for 5 parameters'])
disp([num2str(...
100* wmean(reshape((squeeze(mNRMSE(1:end-2,1,[5],:)) - squeeze(mNRMSE(1:end-2,2,[5],:)))...
    ./ squeeze(mNRMSE(1:end-2,1,[5],:)),1,[]), repmat(w,4,1)')...
, 3) ' +/- ' num2str(...
100* var(reshape((squeeze(mNRMSE(1:end-2,1,[5],:)) - squeeze(mNRMSE(1:end-2,2,[5],:)))...
    ./ squeeze(mNRMSE(1:end-2,1,[5],:)),1,[]), repmat(w,4,1)').^.5...
, 3) ' % for 7 parameters'])


disp(' ')
disp(['For the DBM method the difference between the biggest and the smallest'])
disp(['dictionaries (7 parameters) is ' ...
    num2str(100* wmean(reshape((mRMSE(1:end-2,1,5,1) - mRMSE(1:end-2,1,5,end)) ./ mRMSE(1:end-2,1,5,1),1,[]), w'),3) ...
    ' +/- ' num2str(100* var(reshape((mRMSE(1:end-2,1,5,1) - mRMSE(1:end-2,1,5,end)) ./ mRMSE(1:end-2,1,5,1),1,[]), w').^.5, 3)])
disp(['For the DBL method the difference between the biggest and the smallest'])
disp(['dictionaries (7 parameters) is ' ...
    num2str(100* wmean(reshape((mRMSE(1:end-2,2,5,1) - mRMSE(1:end-2,2,5,end)) ./ mRMSE(1:end-2,2,5,1),1,[]), w'),3) ...
    ' +/- ' num2str(100* var(reshape((mRMSE(1:end-2,2,5,1) - mRMSE(1:end-2,2,5,end)) ./ mRMSE(1:end-2,2,5,1),1,[]), w').^.5, 3)])
disp(['For the DBL method the difference between the 2 highest dictionaries'])
disp(['(7 parameters) is ' ...
    num2str(100* wmean(reshape((mRMSE(1:end-2,2,5,end-1) - mRMSE(1:end-2,2,5,end)) ./ mRMSE(1:end-2,2,5,end-1),1,[]), w'),3) ...
    ' +/- ' num2str(100* var(reshape((mRMSE(1:end-2,2,5,end-1) - mRMSE(1:end-2,2,5,end)) ./ mRMSE(1:end-2,2,5,end-1),1,[]), w').^.5, 3)])


disp(' ')
disp('Altogether, compared to the DBM method, the DBL method reduces the RMSE by:')
disp([num2str(...
100* wmean((squeeze(mNRMSE(1:end-2,1,[3],end)) - squeeze(mNRMSE(1:end-2,2,[3],1)))...
    ./ squeeze(mNRMSE(1:end-2,1,[3],end)), w)...
, 3) ' +/- ' num2str(...
100* var((squeeze(mNRMSE(1:end-2,1,[3],end)) - squeeze(mNRMSE(1:end-2,2,[3],1)))...
    ./ squeeze(mNRMSE(1:end-2,1,[3],end)),w).^.5...
, 3) ' % for 5 parameters'])

disp([num2str(...
100* wmean((squeeze(mNRMSE(1:end-2,1,[5],end)) - squeeze(mNRMSE(1:end-2,2,[5],1)))...
    ./ squeeze(mNRMSE(1:end-2,1,[5],end)), w)...
, 3) ' +/- ' num2str(...
100* var((squeeze(mNRMSE(1:end-2,1,[5],end)) - squeeze(mNRMSE(1:end-2,2,[5],1)))...
    ./ squeeze(mNRMSE(1:end-2,1,[5],end)),w).^.5...
, 3) ' % for 7 parameters'])