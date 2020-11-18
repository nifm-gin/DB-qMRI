
%% Description
%
% We investigate the impact of three parameter sampling strategies
% (regular, random and quasi-random) on the accuracy of parameter
% estimation, using toys signals. We considered different parameter space
% dimensions and dictionary sizes. For each dictionary dimension, test toy
% signals are produced based on a random sampling of the parameter space.
%
% Fabien Boux - 01/2020

Init
disp(['Running experiment ' mfilename '.m'])


%% Setting

% Execution settings
verbose = 2; %0, 1 or 2 for more details
backup  = 1;

% Signal settings
int   	= [.01 1];
p       = [.01 .01];
while min(abs(pdist(p',@(x,y) x-y))) < .05, p = 0.1 + 0.9*rand(1,10); end

param   = [3 3 3 ...
           5 5 5 ...
           7 7];
signals = [512    1728   4096 ... %8^3 12^3 16^3
           1024   3125   7776 ... %4^5 5^5 6^5
           2187  16384];          %3^7 4^7 

% Experiment settings
nb_test_signals = 1000;
nb_repetition = 500;
snr_test = inf;

% Method settings
methods = {'DBM', 'DB-DL'}; %{'DBM', 'DB-SL'};
sampling_strategies = {'Grid', 'Random', 'qRandom'};
Model.snrtrain = inf;
Model.Exec = 'cpu';


%% Creating data

% Expe rand ID
Id      = randi(2^16,1);
% You can specify a value to use an experiment already saved in the temp folder


%% Processing 

if verbose >= 1, disp('Nber of parameters - Dictionary size'); end

for exp = 1:length(signals) % for all experiments
    
    nb_param = param(exp);
    nb_train_signals = signals(exp);
        
    if verbose >= 1, disp([num2str(nb_param) ' - ' num2str(nb_train_signals)]); end
    
    if ~exist(['temp/ParameterSpaceSampling/' num2str(Id)],'dir')
        mkdir(['temp/ParameterSpaceSampling/' num2str(Id)])
    end
    
    if ~exist(['temp/ParameterSpaceSampling/' num2str(Id) '/' num2str(nb_param) '-' num2str(nb_train_signals) '.mat'], 'file')
        
        mRMSE_DBM   = nan(length(sampling_strategies), nb_repetition);
        mMAE_DBM    = mRMSE_DBM; mNRMSE_DBM = mRMSE_DBM; mNMAE_DBM = mRMSE_DBM;
        mRMSE_DBSL  = mRMSE_DBM; mMAE_DBSL = mRMSE_DBM; mNRMSE_DBSL = mRMSE_DBM; mNMAE_DBSL = mRMSE_DBM;
        mRMSE_DBDL  = mRMSE_DBM; mMAE_DBDL = mRMSE_DBM; mNRMSE_DBDL = mRMSE_DBM; mNMAE_DBDL = mRMSE_DBM;
        sizes       = nan(length(sampling_strategies), nb_repetition);

        for rep = 1:nb_repetition

            if verbose == 2, disp([num2str(rep) '/' num2str(nb_repetition)]); end

            % Generate test data
            [Xtest,  Ytest]  = GenerateScalableSignals(p(1:nb_param), int, nb_test_signals, 'Random');
            Xtest   = AddNoise(Xtest, snr_test);

            parfor s = 1:length(sampling_strategies)

                [X, Y]  = GenerateScalableSignals(p(1:nb_param), int, nb_train_signals, sampling_strategies{s});
                sizes(s,rep) = size(Y,1);

                % Prepare dico
                Dico = FormatDico(abs(X), Y);

                % Compute estimates
                if any(contains(methods,'DBM'))
                    tic;
                    Estim       = AnalyzeMRImages(Xtest, Dico, 'DBM', [], Ytest(:,1:size(Y,2)));
                    t_DBM(s,rep)        = toc;
                    
                    mRMSE_DBM(s,rep)	= mean(Estim.GridSearch.Errors.Rmse);
                    mMAE_DBM(s,rep)     = mean(Estim.GridSearch.Errors.Mae);
                    mNRMSE_DBM(s,rep)   = mean(Estim.GridSearch.Errors.Nrmse);
                    mNMAE_DBM(s,rep)    = mean(Estim.GridSearch.Errors.Nmae);
                end

                if any(contains(methods,'DB-SL'))
                    
                    tic;
                    Estim       = AnalyzeMRImages(Xtest, Dico, 'DB-SL', Model, Ytest(:,1:size(Y,2)));
                    t_DBL(s,rep)        = toc;
                    
                    mRMSE_DBSL(s,rep)   = mean(Estim.Regression.Errors.Rmse);
                    mMAE_DBSL(s,rep)    = mean(Estim.Regression.Errors.Mae);
                    mNRMSE_DBSL(s,rep)  = mean(Estim.Regression.Errors.Nrmse);
                    mNMAE_DBSL(s,rep)   = mean(Estim.Regression.Errors.Nmae);
                end
                
                if any(contains(methods,'DB-DL'))
                    
                    tic;
                    Estim       = AnalyzeMRImages(Xtest, Dico, 'DB-DL', Model, Ytest(:,1:size(Y,2)));
                    t_DBDL(s,rep)        = toc;
                    
                    mRMSE_DBDL(s,rep)    = mean(Estim.Regression.Errors.Rmse);
                    mMAE_DBDL(s,rep)     = mean(Estim.Regression.Errors.Mae);
                    mNRMSE_DBDL(s,rep)   = mean(Estim.Regression.Errors.Nrmse);
                    mNMAE_DBDL(s,rep)    = mean(Estim.Regression.Errors.Nmae);
                end
            end   
        end

        % Save intermediate results
        clear *tmp* Dico X* Y* Estim
        save(['temp/ParameterSpaceSampling/' num2str(Id) '/' num2str(nb_param) '-' num2str(nb_train_signals) '.mat'])

        % % Display histogram of errors
        % fig = figure;
        % 
        % subplot(121)
        % hist(mRMSE',0:0.0005:0.08)
        % xlabel('Average RMSE (s) (seconds)'); ylabel('Histogram')
        % legend({'Regular grid', 'Random uniform sampling', 'Quasi random sampling'})
        % xlim([.005 .08]) 
        % 
        % subplot(122)
        % hist(mMAE',0:0.0005:0.08)
        % xlabel('Mean MAE (seconds)'); ylabel('Histogram')
        % legend({'Regular grid', 'Random uniform sampling', 'Quasi random sampling'})
        % xlim([.005 .08])
    end
end

   
%% Displaying

DBM_rmse = struct([]);
DBL_rmse = struct([]);

for i = 1:length(signals)
    filename = ['temp/ParameterSpaceSampling/' num2str(Id) '/' num2str(param(i)) '-' num2str(signals(i)) '.mat'];
    
    if exist(filename,'file')
        load(filename,'mRMSE_DBM','mRMSE_DBSL','mRMSE_DBDL','nb_param','nb_train_signals');
        DBM_rmse{i} = mRMSE_DBM;
        DBL_rmse{i} = mRMSE_DBDL;
    else
        DBM_rmse{i} = [];
        DBL_rmse{i} = [];
    end
    
    title_nb_param(i) = nb_param;
    title_nb_signals(i) = nb_train_signals;
end


% main figure
titles = {'(a)', '(b)', '(c)'};

fig = figure;
c = 0;
for i = [1 4 7]
    c = c + 1;
    
    dat     = DBL_rmse{i}';
    h(c)    = subplot(1,3,c);
    
    boxplot(dat,'symbol','')
    if c == 1
        ylabel('Average RMSE (s)')
    else
        ylabel(' ')
    end
    xtickangle(45)
    set(gca, 'fontsize',15, 'XTickLabels', {'Regular','Rand','qRand'})
    
    title([titles{c} ' ' num2str(title_nb_param(i)) ' parameters'])
    
%     m = mean(dat);
%     disp(['qRand / Grid: ' num2str(1 - m(3) / m(1))])
%     disp(['qRand / Rand: ' num2str(1 - m(3) / m(2))])
end

linkaxes(h,'y')


%% supplementary figure
titles  = {'(a)','(b)','(c)','(d)','(e)','(f)','(g)','(h)'};
alpha   = 0.001;

fig_supp = figure;

for i = 1:length(DBL_rmse)
    
    dat     = [DBM_rmse{i}; DBL_rmse{i}]';
    h(i)    = subplot(3,3,i);
    
    boxplot(dat,'symbol','','OutlierSize',1)
    if i == 1 || i == 4 || i == 7
        ylabel([num2str(title_nb_param(i)) ' parameters \newline Average RMSE (s)'])
    else
        ylabel(' ')
    end
    xtickangle(45)
    set(gca,'linew',1.5, 'fontsize',15, 'XTickLabels',{'Grid','Rand','qRand'})
    
    title(titles{i})
    set(findobj(gca,'type','line'),'linew',1.5)
    xtickangle(45)
    pause(.2)
    
    
    xt  = get(gca, 'XTick');
    hold on
    m1  = max(dat(:,1));
    m2  = max(dat(:,2));
    m3  = max(dat(:,3));
    count = 0;
    boxplot(dat,'symbol','','OutlierSize',1)
    if i == 1 || i == 4 || i == 7
        ylabel([num2str(title_nb_param(i)) ' parameters \newline Average RMSE (s)'])
    else
        ylabel(' ')
    end
    xtickangle(45)
    set(gca,'linew',1.5, 'fontsize',15, 'XTickLabels',{'Regular','Rand','QRand'})
    
    title([titles{c} ' ' num2str(title_nb_signals(i)) ' signals'])
    set(findobj(gca,'type','line'),'linew',1.5)
    xtickangle(45)
    pause(.2)
    
    if verbose >= 2
        disp(['KS test: ' num2str(kstest(dat(:,1),'Alpha',alpha)) ' - ' num2str(kstest(dat(:,2),'Alpha',alpha)) ' - ' num2str(kstest(dat(:,3),'Alpha',alpha))])
        disp(['T test: ' num2str(ttest2(dat(:,1),dat(:,2),'Alpha',alpha,'Vartype','unequal')) ' - ' num2str(ttest2(dat(:,1),dat(:,3),'Alpha',alpha,'Vartype','unequal')) ' - ' num2str(ttest2(dat(:,2),dat(:,3),'Alpha',alpha,'Vartype','unequal'))])
    end
    
    
%     if kstest(dat(:,2),'Alpha',alpha) == 1 && ...
%             kstest(dat(:,3),'Alpha',alpha) == 1 && ...
%             ttest2(dat(:,2),dat(:,3),'Alpha',alpha,'Vartype','unequal') == 1
%         plot(xt([2 3]), 1.05*[max([m2 m3]) max([m2 m3])], '-k',  mean(xt([2 3])), 0.01+1.05*max([m2 m3]), '*k', 'linewidth',1.5, 'markersize',10,'LineWidth', 2)
%         count = count + 1;
%     end
%     if kstest(dat(:,1),'Alpha',alpha) == 1 && ...
%             kstest(dat(:,2),'Alpha',alpha) == 1 && ...
%             ttest2(dat(:,1),dat(:,2),'Alpha',alpha,'Vartype','unequal') == 1
%         plot(xt([1 2]), 1.05*[max([m2 m1]) max([m2 m1])], '-k',  mean(xt([1 2])), 0.01+1.05*max([m2 m1]), '*k', 'linewidth',1.5, 'markersize',10,'LineWidth', 2)
%         count = count + 1;
%     end
%     count = count + .01;
%     if kstest(dat(:,1),'Alpha',alpha) == 1 && ...
%             kstest(dat(:,3),'Alpha',alpha) == 1 && ...
%             ttest2(dat(:,1),dat(:,3),'Alpha',alpha,'Vartype','unequal') == 1
%         if count == 0 
%             plot(xt([1 3]), 1.05*[max([m1 m3]) max([m1 m3])], '-k',  mean(xt([1 3])), 0.01+1.05*max([m1 m3]), '*k', 'linewidth',1.5, 'markersize',10,'LineWidth', 2)
%         else
%             plot(xt([1 3]), 1.3*[max([m1 m3]) max([m1 m3])], '-k',  mean(xt([1 3])), 0.01+1.3*max([m1 m3]), '*k', 'linewidth',1.5, 'markersize',10,'LineWidth', 2)
%         end
%     end
%     hold off
    
%     m = mean(dat);
%     disp(['qRand / Grid: ' num2str(1 - m(6) / m(4))])
%     disp(['qRand / Rand: ' num2str(1 - m(6) / m(5))])
end

linkaxes(h(1:3), 'y')
linkaxes(h(4:6), 'y')
linkaxes(h(7:8), 'y')


%% Sampling strategies illustration

grid_bvf        = 0.05:0.05:0.95;
grid_vsi        = 0.05:0.05:0.95;
[Xgrid, Ygrid]  = meshgrid(grid_bvf, grid_vsi);

nb_signal       = size(Xgrid,1) * size(Xgrid,2);

rand_bvf        = rand(1, nb_signal);
rand_vsi        = rand(1, nb_signal);

vect            = net(scramble(sobolset(2),'MatousekAffineOwen'),nb_signal);
qrand_bvf       = vect(:,1);
qrand_vsi       = vect(:,2);

markersize = 10;
xl  = [0 1];
yl  = [0 1];

fig_illustration = figure;

clear h
h(1) = subplot(131); plot(reshape(Xgrid,1,[]), reshape(Ygrid,1,[]), '.', 'LineWidth', 2, 'MarkerSize', markersize)
xlabel('Second parameter'); ylabel('First parameter')
title('(a)') %title('Regular sampling')
xlim(xl); ylim(yl)

h(2) = subplot(132); plot(rand_bvf, rand_vsi, '.', 'LineWidth', 2, 'MarkerSize', markersize)
xlabel('Second parameter'); ylabel('First parameter')
title('(b)') %title('Random sampling')
xlim(xl); ylim(yl)

h(3) = subplot(133); plot(qrand_bvf, qrand_vsi, '.', 'LineWidth', 2, 'MarkerSize', markersize)
xlabel('Second parameter'); ylabel('First parameter')
title('(c)') %title('quasi-Random sampling')
xlim(xl); ylim(yl)

set(h, 'fontsize', 18, 'XTickLabel',[], 'YTickLabel',[])


%% Exporting figures

if backup == 1
    savefig(fig, 'figures/ParameterSpaceSampling')
    savefig(fig_supp, 'figures/ParameterSpaceSampling-supp')
    savefig(fig_illustration, 'figures/ParameterSpaceSampling-illustration')
end

