
%% Description
%
% We investigate the impact of three parameter sampling strategies
% (regular, random and quasi-random) on the accuracy of parameter
% estimation, using toys signals. We considered different parameter space
% dimensions and dictionary sizes. For each dictionary dimension, test toy
% signals are produced based on a random sampling of the parameter space.
%
% Fabien Boux - 01/2020


%% Setting

% Execution settings
verbose = 1; %0, 1 or 2 for more details
backup  = 1;

% Signal settings
int   	= [0 1];
p       = [.01 .01];
while min(abs(pdist(p',@(x,y) x-y))) < .05, p = 0.1 + 0.9*rand(1,10); end

param   = [3 3 3 ...
           5 5 5 ...
           7 7];
signals = [216  1000 4096 ... %6^3 10^3 16^3
           1024 3125 7776 ... %4^5 5^5 6^5
           2187 16384];       %3^7 4^7 

% Experiment settings

nb_test_signals = 1000;
nb_repetition = 1000;
snr  	= inf;

% Method settings
methods = {'ClassicMRF', 'RegressionMRF'};
Parameters.K = 50;
Parameters.cstr.Sigma  = 'd*';
Parameters.cstr.Gammat = ''; 
Parameters.cstr.Gammaw = '';
Parameters.Lw = 0;
snr_train = inf;


%% Creating data

% Adding to path
addpath(genpath('functions'))
addpath(genpath('tools'))

% Expe rand ID
Id      = randi(2^16,1);
% You can specify a value to replay an experiment saved in the temp folder


%% Processing 

if verbose >= 1, disp('Nber of parameters - Dictionary size'); end

for exp = 1:length(signals) % for all experiments
    
    nb_param = param(exp);
    nb_train_signals = signals(exp);
        
    if verbose >= 1, disp([num2str(nb_param) ' - ' num2str(nb_train_signals)]); end
    
    mRMSE_DBM = nan(3, nb_repetition);
    mMAE_DBM = mRMSE_DBM; mNRMSE_DBM = mRMSE_DBM; mNMAE_DBM = mRMSE_DBM;
    mMAE_DBL = mRMSE_DBM; mNRMSE_DBL = mRMSE_DBM; mNMAE_DBL = mRMSE_DBM;

    parfor rep = 1:nb_repetition

        if verbose == 2, disp([num2str(rep) '/' num2str(nb_repetition)]); end

        % Generate test data
        Ytest   = int(1) + (int(2) - int(1)) * rand(nb_test_signals,nb_param);
        Xtest   = [];
        for sim = 1:size(Ytest,1)
            Xtest(sim,:) = toyMRsignal(Ytest(sim,:),p(1:nb_param));
        end
        Xtest   = AddNoise(Xtest, snr);
        
        for s = 1:3 %for the 3 sampling strategies

            X = []; Y = [];
            switch s 
                % Generate dico grid
                case 1 
                    step    = (int(2)-int(1)) / (nb_train_signals^(1/nb_param));
                    v       = int(1)+step/2:step:int(2)-step/2;
                    Y       = arrangement(v,nb_param);
                    for sim = 1:size(Y,1)
                        X(sim,:) = toyMRsignal(Y(sim,:),p(1:nb_param));
                    end
                    if nb_train_signals ~= length(Y), warning('Not enought signals in grid'); end

                % Generate random uniform dico    
                case 2 
                    Y       = int(1) + (int(2) - int(1)) * rand(nb_train_signals,nb_param);
                    for sim = 1:size(Y,1)
                        X(sim,:) = toyMRsignal(Y(sim,:),p(1:nb_param));
                    end

                % Generate quasi-random dico
                case 3 
                    Y       = int(1) + (int(2) - int(1)) * net(scramble(sobolset(nb_param),'MatousekAffineOwen'),nb_train_signals);
                    for sim = 1:size(Y,1)
                        X(sim,:) = toyMRsignal(Y(sim,:),p(1:nb_param));
                    end
            end

            % Prepare dico
            Dico = [];
            Dico{1}.MRSignals = abs(X); 
            Dico{1}.Parameters.Par = Y;

            % Compute estimates
            if any(contains(methods,'ClassicMRF'))
                tic;
                Estim 	= AnalyzeMRImages(Xtest,Dico,'ClassicMRF',[],Ytest(:,1:size(Dico{1}.Parameters.Par,2)));

                t_DBM(s,rep)        = toc;
                mRMSE_DBM(s,rep)	= mean(Estim.GridSearch.Errors.Rmse);
                mMAE_DBM(s,rep)     = mean(Estim.GridSearch.Errors.Mae);
                mNRMSE_DBM(s,rep)   = mean(Estim.GridSearch.Errors.Nrmse);
                mNMAE_DBM(s,rep)    = mean(Estim.GridSearch.Errors.Nmae);
            end
            
            if any(contains(methods,'RegressionMRF'))
                Dico{1}.MRSignals = AddNoise(Dico{1}.MRSignals, snr_train);

                tic;
                Estim 	= AnalyzeMRImages(Xtest,Dico,'RegressionMRF',Parameters,Ytest(:,1:size(Dico{1}.Parameters.Par,2)));

                t_DBL(s,rep)        =  toc;
                mRMSE_DBL(s,rep)    = mean(Estim.Regression.Errors.Rmse);
                mMAE_DBL(s,rep)     = mean(Estim.Regression.Errors.Mae);
                mNRMSE_DBL(s,rep)   = mean(Estim.Regression.Errors.Nrmse);
                mNMAE_DBL(s,rep)    = mean(Estim.Regression.Errors.Nmae);
            end
        end   
    end

    % Save results
    clear *tmp* Dico X* Y* Estim
    if ~exist(['temp/ParameterSpaceSampling/' num2str(Id)],'dir')
        mkdir(['temp/ParameterSpaceSampling/' num2str(Id)])
    end
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

   
%% Displaying

DBM_rmse = struct([]);
DBL_rmse = struct([]);

for i = 1:length(signals)
    filename = ['temp/ParameterSpaceSampling/' num2str(Id) '/' num2str(param(exp)) '-' num2str(signals(exp)) '.mat'];
    
    if exist(filename,'file')
        load(filename,'mRMSE_DBM','mRMSE_DBL');
        DBM_rmse{i} = mRMSE_DBM;
        DBL_rmse{i} = mRMSE_DBL;
    else
        DBM_rmse{i} = [];
        DBL_rmse{i} = [];
    end
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
    set(gca, 'fontsize',15, 'XTickLabels',{'Grid','Rand','QRand'})
    
    title(titles{c})
end


%% supplementary figure
titles  = {'(a)','(b)','(c)','(d)','(e)','(f)','(g)','(h)'};
alpha   = 0.001;

fig_supp = figure;

for i = 1:length(DBL_rmse)
    
    dat     = [DBM_rmse{i}; DBL_rmse{i}]';
    h(i)    = subplot(3,3,i);
    
    boxplot(dat,'symbol','','OutlierSize',1)
    if i == 1 || i == 4 || i == 7
        ylabel('Average RMSE (s)')
    else
        ylabel(' ')
    end
    xtickangle(45)
    set(gca,'linew',1.5, 'fontsize',15, 'XTickLabels',{'Grid','Rand','QRand'})
    
    title(titles{i})
    set(findobj(gca,'type','line'),'linew',1.5)
    xtickangle(45)
    pause(.2)
    
    if verbose >= 1
%         disp(['KS test: ' num2str(kstest(dat(:,1),'Alpha',alpha)) ' - ' num2str(kstest(dat(:,2),'Alpha',.0001)) ' - ' num2str(kstest(dat(:,3),'Alpha',alpha))])
%         disp(['T test: ' num2str(ttest2(dat(:,1),dat(:,2),'Alpha',alpha,'Vartype','unequal')) ' - ' num2str(ttest2(dat(:,1),dat(:,3),'Alpha',alpha,'Vartype','unequal')) ' - ' num2str(ttest2(dat(:,2),dat(:,3),'Alpha',alpha,'Vartype','unequal'))])
    end
    
    xt  = get(gca, 'XTick');
    hold on
    m1  = max(dat(:,1));
    m2  = max(dat(:,2));
    m3  = max(dat(:,3));
    count = 0;
    boxplot(dat,'symbol','','OutlierSize',1)
    if i == 1 || i == 4 || i == 7
        ylabel('Average RMSE (s)')
    else
        ylabel(' ')
    end
    xtickangle(45)
    set(gca,'linew',1.5, 'fontsize',15, 'XTickLabels',{'Grid','Rand','QRand'})
    
    title(titles{i})
    set(findobj(gca,'type','line'),'linew',1.5)
    xtickangle(45)
    pause(.2)
    
    if verbose >= 1
%         disp(['KS test: ' num2str(kstest(dat(:,1),'Alpha',alpha)) ' - ' num2str(kstest(dat(:,2),'Alpha',.0001)) ' - ' num2str(kstest(dat(:,3),'Alpha',alpha))])
%         disp(['T test: ' num2str(ttest2(dat(:,1),dat(:,2),'Alpha',alpha,'Vartype','unequal')) ' - ' num2str(ttest2(dat(:,1),dat(:,3),'Alpha',alpha,'Vartype','unequal')) ' - ' num2str(ttest2(dat(:,2),dat(:,3),'Alpha',alpha,'Vartype','unequal'))])
    end
    
%     xt  = get(gca, 'XTick');
%     hold on
%     m1 = max(dat(:,1));
%     m2 = max(dat(:,2));
%     m3 = max(dat(:,3));
%     count = 0;
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
    
    if i == 1
        ylim([0 0.2])
    elseif i == 2 || i == 3
        ylim([0.01 0.11])
    elseif i == 4 || i == 5 || i == 6
        ylim([0.03 0.16])
    elseif i > 6
        ylim([0.05 0.2])
    end
end

linkaxes(h(1), 'y')
ylim(h(1),[0.0 0.1])
linkaxes(h(2:3), 'y')
ylim(h(2:3),[0.01 0.11])

linkaxes(h(4:6), 'y')
ylim(h(4:6),[0.0 0.1])

linkaxes(h(7:8), 'y')
ylim(h(7:8),[0.0 0.1])


%% Sampling strategies illustration

grid_bvf        = 0.05:0.05:0.95;
grid_vsi        = 0.05:0.05:0.95;
[Xgrid, Ygrid]  = meshgrid(grid_bvf, grid_vsi);

nb_signal       = size(Xgrid,1) * size(Xgrid,2);

rand_bvf        = rand(1, nb_signal);
rand_vsi        = rand(1, nb_signal);

% vect            = net(sobolset(2),nb_signal);
vect            = net(scramble(sobolset(2),'MatousekAffineOwen'),nb_signal);
qrand_bvf       = vect(:,1);
qrand_vsi       = vect(:,2);

markersize = 10;
xl = [0 1];
yl = [0 1];

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

