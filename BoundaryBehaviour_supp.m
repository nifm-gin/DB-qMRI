
%% Description
%
% The DBL method estimates the parameters using a continuous function that
% is not restricted to the dictionary parameter space. We therefore 
% investigated the behaviors of the DBL and DBM methods outside the
% boundaries of this space. Two 2-dimensional parameter subspaces was
% defined: one to sample the parameters of the dictionary, and a larger one
% for the test.
% Compared to the BoundaryBehaviour.m script, we investigate the impact of
% adding few signals on the estimates.
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
nb_test_signals =  100;
nb_train_signals = 1000;
N     	= 198;
lw      = 10;
nb_int  = 2; % Can be 2 or 3 (default 2)
int_all = {{[.20 .30], [.30 .60]},...
           {[.70 .80], [.10 .25]},...
           {[.40 .65], [.30 .45]},...
           };
nb      = 3; % number of added signals

% Regression settings
Parameters = [];
Parameters.K = 50;
Parameters.cstr.Sigma  = 'd*';
Parameters.cstr.Gammat = ''; 
Parameters.cstr.Gammaw = '';
Parameters.Lw = 0;
snr_train  	= 100;


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

% Add some values out of the training space
i       = i+1;
% newY    = int_all{i}{1}(1) + (int_all{i}{1}(2) - int_all{i}{1}(1)) * net(scramble(sobolset(nb_param),'MatousekAffineOwen'),nb);
newY = [0.5689    0.5159
    0.8333    0.7426
    0.4917    0.7606];
newX    = [];
for sim = 1:size(newY,1), newX(sim,:) = toyMRsignal(newY(sim,:), p); end
Ytrain      = [Ytrain; newY];
Xtrain      = [Xtrain; newX];
Ytrain_grid = [Ytrain_grid; newY];
Xtrain_grid = [Xtrain_grid; newX];


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

    parfor s2 = floor(lw/2)+1:length(intt_)-floor(lw/2)

        % Test data
        subint1 = intt_([s1-floor(lw/2) s1+floor(lw/2)]);
        subint2 = intt_([s2-floor(lw/2) s2+floor(lw/2)]);
        inter2(s2) = intt_(s2);

        Ytest = [];
        tmp         = net(scramble(sobolset(nb_param),'MatousekAffineOwen'),nb_test_signals);
        Ytest(:,1)  = subint1(1) + (subint1(2) - subint1(1)) * tmp(:,1);
        Ytest(:,2)  = subint2(1) + (subint2(2) - subint2(1)) * tmp(:,2);

        Xtest = [];
        for sim = 1:size(Ytest,1)
            Xtest(sim,:) = toyMRsignal(Ytest(sim,:), p);
        end
        Xtest_  = AddNoise(Xtest, snr);

        % Compute estimates
        Estim   = AnalyzeMRImages(Xtest_, [], 'DBL', Parameters); 
        Ygllim  = Estim.Regression.Y(:,1:nb_param);
        [Rmse_gllim(s1,s2,:),~, Mae_gllim(s1,s2,:)] = EvaluateEstimation(Ytest, Ygllim);

        Ygrid   = EstimateParametersFromGrid(Xtest_,Xtrain_grid,Ytrain_grid);
        [Rmse_grid(s1,s2,:),~, Mae_grid(s1,s2,:)] = EvaluateEstimation(Ytest, Ygrid);
    end
end
    

%% Saving 

if backup == 1
    clear tmp* X* Y* Dic
    save(['temp/' 'BoundaryBehaviour-supp'])
end


%% Display

fig = figure;

v       = 1:24:length(intt_);
ttt     = split(num2str(intt_(v),1), ' ');
ttt     = ttt(~cellfun('isempty',ttt));

err     = Rmse_grid;
bounds  = [0 .2];

h(nb_int-1) = subplot(2,length(int_all)-1,nb_int-1);

imagesc(inter1,inter2,mean(err,3), bounds)
hold on
for i = 1:length(int)
    line([int{i}{2}(1) int{i}{2}(2)], [int{i}{1}(1) int{i}{1}(1)], 'linestyle', '--', 'color','w',  'linewidth',3)
    line([int{i}{2}(1) int{i}{2}(2)], [int{i}{1}(2) int{i}{1}(2)], 'linestyle', '--', 'color','w',  'linewidth',3)
    line([int{i}{2}(1) int{i}{2}(1)], [int{i}{1}(1) int{i}{1}(2)], 'linestyle', '--', 'color','w',  'linewidth',3)
    line([int{i}{2}(2) int{i}{2}(2)], [int{i}{1}(1) int{i}{1}(2)], 'linestyle', '--', 'color','w',  'linewidth',3)
end
plot(newY(:,2),newY(:,1), 'wx', 'markersize', 12)    
colormap(mycmap()); colorbar
xlabel('Second parameter'); ylabel('First parameter')
xlim([intt_(floor(lw/2)+1) intt_(length(intt_)-floor(lw/2))])
ylim([intt_(floor(lw/2)+1) intt_(length(intt_)-floor(lw/2))])

if nb_int == 2
    title('(a)')
elseif nb_int == 3
    title('(b)')
end
set(gca, 'fontsize', 18)
set(gca,'DataAspectRatio',[10 10 10])
set(gca,'YDir','normal')

h(length(int_all)-1+nb_int) = subplot(2,length(int_all)-1,length(int_all)-1+nb_int-1);
err     = Rmse_gllim;
imagesc(inter1,inter2,mean(err,3), bounds)
hold on
for i = 1:length(int)
    line([int{i}{2}(1) int{i}{2}(2)], [int{i}{1}(1) int{i}{1}(1)], 'linestyle', '--', 'color','w',  'linewidth',3)
    line([int{i}{2}(1) int{i}{2}(2)], [int{i}{1}(2) int{i}{1}(2)], 'linestyle', '--', 'color','w',  'linewidth',3)
    line([int{i}{2}(1) int{i}{2}(1)], [int{i}{1}(1) int{i}{1}(2)], 'linestyle', '--', 'color','w',  'linewidth',3)
    line([int{i}{2}(2) int{i}{2}(2)], [int{i}{1}(1) int{i}{1}(2)], 'linestyle', '--', 'color','w',  'linewidth',3)
end
plot(newY(:,2),newY(:,1), 'wx', 'markersize', 12)    
colormap(mycmap()); colorbar
xlabel('Second parameter'); ylabel('First parameter')
xlim([intt_(floor(lw/2)+1) intt_(length(intt_)-floor(lw/2))])
ylim([intt_(floor(lw/2)+1) intt_(length(intt_)-floor(lw/2))])

if nb_int == 2
    title('(c)')
elseif nb_int == 3
    title('(d)')
end
set(gca, 'fontsize', 18)
set(gca,'DataAspectRatio',[10 10 10])
set(gca,'YDir','normal')


%% Exporting figures

if backup == 1
    savefig(fig, ['figures/' 'BoundaryBehaviour-supp'])
end

