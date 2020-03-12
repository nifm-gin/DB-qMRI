
%% Description
%
% The DBL method estimates the parameters using a continuous function that
% is not restricted to the dictionary parameter space. We therefore 
% investigated the behaviors of the DBL and DBM methods outside the
% boundaries of this space. Two 2-dimensional parameter subspaces was
% defined: one to sample the parameters of the dictionary, and a larger one
% for the test.
%
% Fabien Boux - 01/2020


%% Setting

% Execution settings
verbose = 1;
backup  = 1;

% Signal settings
nb_param = 2;
intt  	= [0.01 1];
p   	= [.5 .8];
snr     = inf;

% Experiment settings
nb_test_signals = 2e6;
nb_train_signals = 5000;
N     	= 198;
lw      = 10;
nb_int  = 2; % Can be 2 or 3 (default 2)

int_all = {{[.20 .30], [.30 .60]},...
           {[.70 .80], [.10 .25]},...
           {[.40 .65], [.30 .45]},...
           };

% Regression settings
Parameters = [];
Parameters.K = 50;
Parameters.cstr.Sigma  = 'd*';
Parameters.cstr.Gammat = ''; 
Parameters.cstr.Gammaw = '';
Parameters.Lw = 0;
snr_train  	= 60;


%% Creating data

% Adding to path
addpath(genpath('functions'))
addpath(genpath('tools'))

% Creating 
int = int_all(1:nb_int);
    
full = 0;
for i = 1:length(int)
    full = full + (int{i}{1}(2) - int{i}{1}(1)) * (int{i}{2}(2) - int{i}{2}(1));
end

% Create training data
clear Xtrain* Ytrain*
for i = 1:length(int)
    proportion = (int{i}{1}(2) - int{i}{1}(1)) * (int{i}{2}(2) - int{i}{2}(1));
    tmp = net(scramble(sobolset(nb_param),'MatousekAffineOwen'),floor(proportion / full * nb_train_signals));
    YtrainS{i}(:,1) = int{i}{1}(1) + (int{i}{1}(2) - int{i}{1}(1)) * tmp(:,1);
    YtrainS{i}(:,2) = int{i}{2}(1) + (int{i}{2}(2) - int{i}{2}(1)) * tmp(:,2);

    ds      = ( (int{i}{1}(2)-int{i}{1}(1)) * (int{i}{2}(2)-int{i}{2}(1)) / ...
                 floor(proportion / full * nb_train_signals) )^.5;

    v1      = int{i}{1}(1)+ds/2:ds:int{i}{1}(2)-ds/2;
    v2      = int{i}{2}(1)+ds/2:ds:int{i}{2}(2)-ds/2;

    [tmpx, tmpy] = meshgrid(v1,v2);
    Ytrain_gS{i} = [tmpx(:) tmpy(:)];

    for sim = 1:size(YtrainS{i},1)
        XtrainS{i}(sim,:) = toyMRsignal(YtrainS{i}(sim,:), p);
    end
    for sim = 1:size(Ytrain_gS{i},1)
        Xtrain_gS{i}(sim,:) = toyMRsignal(Ytrain_gS{i}(sim,:), p);
    end
end
Xtrain = []; Ytrain = [];
Xtrain_grid = []; Ytrain_grid = [];
for i = 1:length(int)
    Ytrain      = [Ytrain; YtrainS{i}];
    Xtrain      = [Xtrain; XtrainS{i}];
    Xtrain_grid = [Xtrain_grid; Xtrain_gS{i}];
    Ytrain_grid = [Ytrain_grid; Ytrain_gS{i}];
end


% Create test signals
Ytest_ = [];
tmp         = rand(nb_test_signals, nb_param);
Ytest_(:,1)  = intt(1) + (intt(2) - intt(1)) * tmp(:,1);
Ytest_(:,2)  = intt(1) + (intt(2) - intt(1)) * tmp(:,2);

Xtest_ = [];
parfor sim = 1:size(Ytest_,1)
    Xtest_(sim,:) = toyMRsignal(Ytest_(sim,:), p);
end
Xtest_  = AddNoise(Xtest_, snr);


%% Processing 

% GLLiM learning
Dic{1}.MRSignals       = AddNoise(Xtrain, snr_train);
Dic{1}.Parameters.Par  = Ytrain;
[~,Parameters]  = AnalyzeMRImages([],Dic,'DBL',Parameters);

% Estimation
intt_   = intt(1):(intt(end) - intt(1))/N:intt(2);

for s1 = floor(lw/2)+1:length(intt_)-floor(lw/2)
    
    if verbose == 1, disp([num2str(s1-floor(lw/2)) '/' num2str(length(floor(lw/2)+1:length(intt_)-floor(lw/2)))]); end

    inter1(s1) = intt_(s1);
    
    subint1 = intt_([s1-floor(lw/2) s1+floor(lw/2)]);
    
    v1      = (subint1(1) <= Ytest_(:,1)) & (Ytest_(:,1) < subint1(2));
    Xtest__ = Xtest_(v1,:);
    Ytest__ = Ytest_(v1,:);
    
    parfor s2 = floor(lw/2)+1:length(intt_)-floor(lw/2)

        % Test data
        
        subint2 = intt_([s2-floor(lw/2) s2+floor(lw/2)]);
        inter2(s2) = intt_(s2);

        v2  	= (subint2(1) <= Ytest__(:,2)) & (Ytest__(:,2) < subint2(2));
        
        Xtest   = Xtest__(v2,:);
        Ytest   = Ytest__(v2,:);
        dens(s1,s2) = sum(v2(:));
        
        % Compute estimates
        Estim   = AnalyzeMRImages(Xtest, [], 'DBL', Parameters); 
        Ygllim  = Estim.Regression.Y(:,1:nb_param);
        [Rmse_gllim(s1,s2,:),~, Mae_gllim(s1,s2,:)] = EvaluateEstimation(Ytest, Ygllim);
        CI(s1,s2,:) = nanmean(squeeze(Estim.Regression.Cov.^.5));
        
        Ygrid   = EstimateParametersFromGrid(Xtest,Xtrain_grid,Ytrain_grid);
        [Rmse_grid(s1,s2,:),~, Mae_grid(s1,s2,:)] = EvaluateEstimation(Ytest, Ygrid);
    end
end
    

%% Saving 

if backup == 1
    clear tmp* X* Y* Dic
    save(['temp/' 'BoundaryBehaviour'])
end


%% Display

fig = figure;

v2       = 1:24:length(intt_);
ttt     = split(num2str(intt_(v2),1), ' ');
ttt     = ttt(~cellfun('isempty',ttt));

err     = Rmse_grid;
bounds = [0 .2];


h(1) = subplot(221);

imagesc(inter1,inter2,mean(err,3), bounds)
hold on
for i = 1:length(int)
    line([int{i}{2}(1) int{i}{2}(2)], [int{i}{1}(1) int{i}{1}(1)], 'linestyle', '--', 'color','w',  'linewidth',3)
    line([int{i}{2}(1) int{i}{2}(2)], [int{i}{1}(2) int{i}{1}(2)], 'linestyle', '--', 'color','w',  'linewidth',3)
    line([int{i}{2}(1) int{i}{2}(1)], [int{i}{1}(1) int{i}{1}(2)], 'linestyle', '--', 'color','w',  'linewidth',3)
    line([int{i}{2}(2) int{i}{2}(2)], [int{i}{1}(1) int{i}{1}(2)], 'linestyle', '--', 'color','w',  'linewidth',3)
end
colormap(mycmap()); colorbar
xlabel('Second parameter'); ylabel('First parameter')
xlim([intt_(floor(lw/2)+1) intt_(length(intt_)-floor(lw/2))])
ylim([intt_(floor(lw/2)+1) intt_(length(intt_)-floor(lw/2))])
title('(a)')

set(gca, 'fontsize', 18)
set(gca,'DataAspectRatio',[10 10 10])
set(gca,'YDir','normal')


h(3) = subplot(223);
err     = Rmse_gllim;
imagesc(inter1,inter2,mean(err,3), bounds)
hold on
for i = 1:length(int)
    line([int{i}{2}(1) int{i}{2}(2)], [int{i}{1}(1) int{i}{1}(1)], 'linestyle', '--', 'color','w',  'linewidth',3)
    line([int{i}{2}(1) int{i}{2}(2)], [int{i}{1}(2) int{i}{1}(2)], 'linestyle', '--', 'color','w',  'linewidth',3)
    line([int{i}{2}(1) int{i}{2}(1)], [int{i}{1}(1) int{i}{1}(2)], 'linestyle', '--', 'color','w',  'linewidth',3)
    line([int{i}{2}(2) int{i}{2}(2)], [int{i}{1}(1) int{i}{1}(2)], 'linestyle', '--', 'color','w',  'linewidth',3)
end
colormap(mycmap()); colorbar
xlabel('Second parameter'); ylabel('First parameter')
xlim([intt_(floor(lw/2)+1) intt_(length(intt_)-floor(lw/2))])
ylim([intt_(floor(lw/2)+1) intt_(length(intt_)-floor(lw/2))])
title('(c)')
    
set(gca, 'fontsize', 18)
set(gca,'DataAspectRatio',[10 10 10])
set(gca,'YDir','normal')


% h(4) = subplot(224);
% 
% imagesc(inter1,inter2,mean(CI,3), bounds/20)
% hold on
% for i = 1:length(int)
%     line([int{i}{2}(1) int{i}{2}(2)], [int{i}{1}(1) int{i}{1}(1)], 'linestyle', '--', 'color','w',  'linewidth',3)
%     line([int{i}{2}(1) int{i}{2}(2)], [int{i}{1}(2) int{i}{1}(2)], 'linestyle', '--', 'color','w',  'linewidth',3)
%     line([int{i}{2}(1) int{i}{2}(1)], [int{i}{1}(1) int{i}{1}(2)], 'linestyle', '--', 'color','w',  'linewidth',3)
%     line([int{i}{2}(2) int{i}{2}(2)], [int{i}{1}(1) int{i}{1}(2)], 'linestyle', '--', 'color','w',  'linewidth',3)
% end
% colormap(mycmap()); colorbar
% xlabel('Second parameter'); ylabel('First parameter')
% xlim([intt_(floor(lw/2)+1) intt_(length(intt_)-floor(lw/2))])
% ylim([intt_(floor(lw/2)+1) intt_(length(intt_)-floor(lw/2))])
% title('(d)')


set(gca, 'fontsize', 18)
set(gca,'DataAspectRatio',[10 10 10])
set(gca,'YDir','normal')


%% Exporting figures

if backup == 1
    savefig(fig, ['figures/' 'BoundaryBehaviour'])
end

