
%% Description
%
% The DBL methods estimate the parameters using a continuous function that
% is not restricted to the dictionary parameter space. We therefore 
% investigated the behaviors of the DBL and DBM methods outside the
% boundaries of this space. Two 2-dimensional parameter subspaces was
% defined: one to sample the parameters of the dictionary, and a larger one
% for the test.
% In addition, we investigate the impact of adding few signals on the
% estimate accuracies.
%
% Fabien Boux - 01/2020

Init
disp(['Running experiment ' mfilename '.m'])


%% Setting

% Execution settings
verbose = 1;
backup  = 1;

% Signal settings
nb_param = 2;
intt  	= [0.01 1];
p   	= [.5 .8];
snr     = 40;

% Experiment settings
nb_test_signals = 2e6;
nb_train_signals = 10000;
N     	= 198;
lw      = 10;
nb_int  = 2; % Can be 2 or 3 (default 2)
int_all = {{[.20 .30], [.30 .60]},...
           {[.70 .80], [.10 .25]},...
           {[.40 .65], [.30 .45]},...
           };


%% Creating data

% Creating 
int = int_all(1:nb_int);
    
full = 0;
for i = 1:length(int)
    full = full + (int{i}{1}(2) - int{i}{1}(1)) * (int{i}{2}(2) - int{i}{2}(1));
end

% Create training data
clear Xtrain* Ytrain*
for i = 1:length(int)
    
    % Design sampling in patches
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
    
    % Simulate
    XtrainS{i} = ScalableSignals(YtrainS{i}, p);
    Xtrain_gS{i} = ScalableSignals(Ytrain_gS{i}, p);
end
Xtrain = []; Ytrain = [];
Xtrain_grid = []; Ytrain_grid = [];
for i = 1:length(int)
    Ytrain      = [Ytrain; YtrainS{i}];
    Xtrain      = [Xtrain; XtrainS{i}];
    Xtrain_grid = [Xtrain_grid; Xtrain_gS{i}];
    Ytrain_grid = [Ytrain_grid; Ytrain_gS{i}];
end

% Add extra dictionary entries out of the training space
i       = i+1;
newY    = [0.5689    0.5159
           0.8333    0.7426
           0.4917    0.7606];
newX    = ScalableSignals(newY, p);

Ytrain_ext      = [Ytrain; newY];
Xtrain_ext      = [Xtrain; newX];
Ytrain_grid_ext = [Ytrain_grid; newY];
Xtrain_grid_ext = [Xtrain_grid; newX];

% Create test signals
[Xtest_, Ytest_]  = GenerateScalableSignals(p, intt, nb_test_signals, 'Random');
Xtest_  = AddNoise(Xtest_, snr);


%% Learnings

% DB-SL/DB-DL learning for initial dico
Dico = FormatDico(Xtrain, Ytrain);
[~,Model]       = AnalyzeMRImages([], Dico, 'DB-SL');
[~,NeuralNet]	= AnalyzeMRImages([], Dico, 'DB-DL');

% DB-SL/DB-DL learning for extended dico (ie with extra entries)
Dico = FormatDico(Xtrain_ext, Ytrain_ext);
[~,Model_ext] 	= AnalyzeMRImages([], Dico, 'DB-SL');
[~,NeuralNet_ext] = AnalyzeMRImages([], Dico, 'DB-DL');


%% Estimation

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
        
        for i = 1:nb_param
            mask(s1,s2,i)   = (( (int{i}{1}(1) ) <= subint1(1) ) &  ( subint1(2) <= int{i}{1}(2) )) & ...
                            (( int{i}{2}(1) <= subint2(1) ) &  ( subint2(2) <= int{i}{2}(2) ));
            mask_out(s1,s2,i)   = (( (int{i}{1}(2) ) <= subint1(1) ) | ( subint1(2) <= int{i}{1}(1) )) | ...
                            (( int{i}{2}(2) <= subint2(1) ) - ( subint2(2) <= int{i}{2}(1) ));
        end
        
        if dens(s1,s2) ~= 0
            % Compute estimates
            %DBM
            Ygrid   = EstimateParametersFromGrid(Xtest,Xtrain_grid, Ytrain_grid);
            [Rmse_grid(s1,s2,:),~, Mae_grid(s1,s2,:)] = EvaluateEstimation(Ytest, Ygrid);
            
            %DB-SL
            Estim   = AnalyzeMRImages(Xtest, [], 'DB-SL', Model, [],[], snr); 
            Ygllim  = Estim.Regression.Y(:,1:nb_param);
            [Rmse_gllim(s1,s2,:),~, Mae_gllim(s1,s2,:)] = EvaluateEstimation(Ytest, Ygllim);
            CI(s1,s2,:) = nanmean(squeeze(Estim.Regression.Cov.^.5));
            %Mahaldist(s1,s2) = nanmean(Estim.Regression.Mahaldist);

            %DB-DL
            Estim   = AnalyzeMRImages(Xtest, [], 'DB-DL', NeuralNet); 
            Ynn     = Estim.Regression.Y(:,1:nb_param);
            [Rmse_nn(s1,s2,:),~, Mae_nn(s1,s2,:)] = EvaluateEstimation(Ytest,Ynn);

            % Compute estimates for extended dico
            %DBM
            Ygrid   = EstimateParametersFromGrid(Xtest,Xtrain_grid_ext, Ytrain_grid_ext);
            [Rmse_grid_ext(s1,s2,:),~, Mae_grid_ext(s1,s2,:)] = EvaluateEstimation(Ytest, Ygrid);
            
            %DB-SL
            Estim   = AnalyzeMRImages(Xtest, [], 'DB-SL', Model_ext, [],[], snr); 
            Ygllim  = Estim.Regression.Y(:,1:nb_param);
            [Rmse_gllim_ext(s1,s2,:),~, Mae_gllim_ext(s1,s2,:)] = EvaluateEstimation(Ytest, Ygllim);
            CI_ext(s1,s2,:) = nanmean(squeeze(Estim.Regression.Cov.^.5));

            %DB-DL
            Estim   = AnalyzeMRImages(Xtest, [], 'DB-DL', NeuralNet_ext); 
            Ynn     = Estim.Regression.Y(:,1:nb_param);
            [Rmse_nn_ext(s1,s2,:),~, Mae_nn_ext(s1,s2,:)] = EvaluateEstimation(Ytest,Ynn);
        end
    end
end
    

%% Saving 

if backup == 1
    clear tmp* X* Y* Dico*
    save(['temp/' 'BoundaryBehavior'])
end


%% Display

fig = figure;

v       = 1:24:length(intt_);
ttt     = split(num2str(intt_(v),1), ' ');
ttt     = ttt(~cellfun('isempty',ttt));

mycmap = mycmap_extended;

bounds  = [0 .2];

h(1) = subplot(341);
err     = Rmse_grid;
imagesc(inter1,inter2,mean(err,3), bounds)
hold on
for i = 1:length(int)
    line([int{i}{2}(1) int{i}{2}(2)], [int{i}{1}(1) int{i}{1}(1)], 'linestyle', '--', 'color','w',  'linewidth',3,'HandleVisibility','off')
    line([int{i}{2}(1) int{i}{2}(2)], [int{i}{1}(2) int{i}{1}(2)], 'linestyle', '--', 'color','w',  'linewidth',3,'HandleVisibility','off')
    line([int{i}{2}(1) int{i}{2}(1)], [int{i}{1}(1) int{i}{1}(2)], 'linestyle', '--', 'color','w',  'linewidth',3,'HandleVisibility','off')
    line([int{i}{2}(2) int{i}{2}(2)], [int{i}{1}(1) int{i}{1}(2)], 'linestyle', '--', 'color','w',  'linewidth',3,'HandleVisibility','off')
end 
colormap(mycmap()); colorbar
xlabel('Second parameter'); ylabel('First parameter')
xlim([intt_(floor(lw/2)+1) intt_(length(intt_)-floor(lw/2))])
ylim([intt_(floor(lw/2)+1) intt_(length(intt_)-floor(lw/2))])
title('(a)')

set(gca, 'fontsize', 18)
set(gca,'DataAspectRatio',[10 10 10])
set(gca,'YDir','normal')



h(3) = subplot(342);
err     = Rmse_nn;
imagesc(inter1,inter2,mean(err,3), bounds)
hold on
for i = 1:length(int)
    line([int{i}{2}(1) int{i}{2}(2)], [int{i}{1}(1) int{i}{1}(1)], 'linestyle', '--', 'color','w',  'linewidth',3,'HandleVisibility','off')
    line([int{i}{2}(1) int{i}{2}(2)], [int{i}{1}(2) int{i}{1}(2)], 'linestyle', '--', 'color','w',  'linewidth',3,'HandleVisibility','off')
    line([int{i}{2}(1) int{i}{2}(1)], [int{i}{1}(1) int{i}{1}(2)], 'linestyle', '--', 'color','w',  'linewidth',3,'HandleVisibility','off')
    line([int{i}{2}(2) int{i}{2}(2)], [int{i}{1}(1) int{i}{1}(2)], 'linestyle', '--', 'color','w',  'linewidth',3,'HandleVisibility','off')
end 
colormap(mycmap()); colorbar
xlabel('Second parameter'); ylabel('First parameter')
xlim([intt_(floor(lw/2)+1) intt_(length(intt_)-floor(lw/2))])
ylim([intt_(floor(lw/2)+1) intt_(length(intt_)-floor(lw/2))])
title('(b)')
    
set(gca, 'fontsize', 18)
set(gca,'DataAspectRatio',[10 10 10])
set(gca,'YDir','normal')



h(3) = subplot(343);
err     = Rmse_gllim;
imagesc(inter1,inter2,mean(err,3), bounds)
hold on
for i = 1:length(int)
    line([int{i}{2}(1) int{i}{2}(2)], [int{i}{1}(1) int{i}{1}(1)], 'linestyle', '--', 'color','w',  'linewidth',3,'HandleVisibility','off')
    line([int{i}{2}(1) int{i}{2}(2)], [int{i}{1}(2) int{i}{1}(2)], 'linestyle', '--', 'color','w',  'linewidth',3,'HandleVisibility','off')
    line([int{i}{2}(1) int{i}{2}(1)], [int{i}{1}(1) int{i}{1}(2)], 'linestyle', '--', 'color','w',  'linewidth',3,'HandleVisibility','off')
    line([int{i}{2}(2) int{i}{2}(2)], [int{i}{1}(1) int{i}{1}(2)], 'linestyle', '--', 'color','w',  'linewidth',3,'HandleVisibility','off')
end
colormap(mycmap()); colorbar
xlabel('Second parameter'); ylabel('First parameter')
xlim([intt_(floor(lw/2)+1) intt_(length(intt_)-floor(lw/2))])
ylim([intt_(floor(lw/2)+1) intt_(length(intt_)-floor(lw/2))])
title('(c)')
    
set(gca, 'fontsize', 18)
set(gca,'DataAspectRatio',[10 10 10])
set(gca,'YDir','normal')





h(5) = subplot(345);
err     = Rmse_grid_ext;
imagesc(inter1,inter2,mean(err,3), bounds)
hold on
for i = 1:length(int)
    line([int{i}{2}(1) int{i}{2}(2)], [int{i}{1}(1) int{i}{1}(1)], 'linestyle', '--', 'color','w',  'linewidth',3,'HandleVisibility','off')
    line([int{i}{2}(1) int{i}{2}(2)], [int{i}{1}(2) int{i}{1}(2)], 'linestyle', '--', 'color','w',  'linewidth',3,'HandleVisibility','off')
    line([int{i}{2}(1) int{i}{2}(1)], [int{i}{1}(1) int{i}{1}(2)], 'linestyle', '--', 'color','w',  'linewidth',3,'HandleVisibility','off')
    line([int{i}{2}(2) int{i}{2}(2)], [int{i}{1}(1) int{i}{1}(2)], 'linestyle', '--', 'color','w',  'linewidth',3,'HandleVisibility','off')
end
plot(newY(:,2),newY(:,1), 'wx', 'markersize', 12,'HandleVisibility','off')
colormap(mycmap()); colorbar
xlabel('Second parameter'); ylabel('First parameter')
xlim([intt_(floor(lw/2)+1) intt_(length(intt_)-floor(lw/2))])
ylim([intt_(floor(lw/2)+1) intt_(length(intt_)-floor(lw/2))])
title('(a)')

set(gca, 'fontsize', 18)
set(gca,'DataAspectRatio',[10 10 10])
set(gca,'YDir','normal')



h(6) = subplot(346);
err     = Rmse_nn_ext;
imagesc(inter1,inter2,mean(err,3), bounds)
hold on
for i = 1:length(int)
    line([int{i}{2}(1) int{i}{2}(2)], [int{i}{1}(1) int{i}{1}(1)], 'linestyle', '--', 'color','w',  'linewidth',3,'HandleVisibility','off')
    line([int{i}{2}(1) int{i}{2}(2)], [int{i}{1}(2) int{i}{1}(2)], 'linestyle', '--', 'color','w',  'linewidth',3,'HandleVisibility','off')
    line([int{i}{2}(1) int{i}{2}(1)], [int{i}{1}(1) int{i}{1}(2)], 'linestyle', '--', 'color','w',  'linewidth',3,'HandleVisibility','off')
    line([int{i}{2}(2) int{i}{2}(2)], [int{i}{1}(1) int{i}{1}(2)], 'linestyle', '--', 'color','w',  'linewidth',3,'HandleVisibility','off')
end
plot(newY(:,2),newY(:,1), 'wx', 'markersize', 12,'HandleVisibility','off')
colormap(mycmap()); colorbar
xlabel('Second parameter'); ylabel('First parameter')
xlim([intt_(floor(lw/2)+1) intt_(length(intt_)-floor(lw/2))])
ylim([intt_(floor(lw/2)+1) intt_(length(intt_)-floor(lw/2))])
title('(b)')
    
set(gca, 'fontsize', 18)
set(gca,'DataAspectRatio',[10 10 10])
set(gca,'YDir','normal')



h(7) = subplot(347);
err     = Rmse_gllim_ext;
imagesc(inter1,inter2,mean(err,3), bounds)
hold on
for i = 1:length(int)
    line([int{i}{2}(1) int{i}{2}(2)], [int{i}{1}(1) int{i}{1}(1)], 'linestyle', '--', 'color','w',  'linewidth',3,'HandleVisibility','off')
    line([int{i}{2}(1) int{i}{2}(2)], [int{i}{1}(2) int{i}{1}(2)], 'linestyle', '--', 'color','w',  'linewidth',3,'HandleVisibility','off')
    line([int{i}{2}(1) int{i}{2}(1)], [int{i}{1}(1) int{i}{1}(2)], 'linestyle', '--', 'color','w',  'linewidth',3,'HandleVisibility','off')
    line([int{i}{2}(2) int{i}{2}(2)], [int{i}{1}(1) int{i}{1}(2)], 'linestyle', '--', 'color','w',  'linewidth',3,'HandleVisibility','off')
end
plot(newY(:,2),newY(:,1), 'wx', 'markersize', 12,'HandleVisibility','off')
colormap(mycmap()); colorbar
xlabel('Second parameter'); ylabel('First parameter')
xlim([intt_(floor(lw/2)+1) intt_(length(intt_)-floor(lw/2))])
ylim([intt_(floor(lw/2)+1) intt_(length(intt_)-floor(lw/2))])
title('(c)')
    
set(gca, 'fontsize', 18)
set(gca,'DataAspectRatio',[10 10 10])
set(gca,'YDir','normal')




h(4) = subplot(344);
err     = CI*10;
imagesc(inter1,inter2,mean(err,3), bounds)
hold on
for i = 1:length(int)
    line([int{i}{2}(1) int{i}{2}(2)], [int{i}{1}(1) int{i}{1}(1)], 'linestyle', '--', 'color','w',  'linewidth',3,'HandleVisibility','off')
    line([int{i}{2}(1) int{i}{2}(2)], [int{i}{1}(2) int{i}{1}(2)], 'linestyle', '--', 'color','w',  'linewidth',3,'HandleVisibility','off')
    line([int{i}{2}(1) int{i}{2}(1)], [int{i}{1}(1) int{i}{1}(2)], 'linestyle', '--', 'color','w',  'linewidth',3,'HandleVisibility','off')
    line([int{i}{2}(2) int{i}{2}(2)], [int{i}{1}(1) int{i}{1}(2)], 'linestyle', '--', 'color','w',  'linewidth',3,'HandleVisibility','off')
end
colormap(mycmap()); colorbar
xlabel('Second parameter'); ylabel('First parameter')
xlim([intt_(floor(lw/2)+1) intt_(length(intt_)-floor(lw/2))])
ylim([intt_(floor(lw/2)+1) intt_(length(intt_)-floor(lw/2))])
title('(supp)')
     
set(gca, 'fontsize', 18)
set(gca,'DataAspectRatio',[10 10 10])
set(gca,'YDir','normal')

h(8) = subplot(348);
err     = CI_ext*10;
imagesc(inter1,inter2,mean(err,3), bounds)
hold on
for i = 1:length(int)
    line([int{i}{2}(1) int{i}{2}(2)], [int{i}{1}(1) int{i}{1}(1)], 'linestyle', '--', 'color','w',  'linewidth',3,'HandleVisibility','off')
    line([int{i}{2}(1) int{i}{2}(2)], [int{i}{1}(2) int{i}{1}(2)], 'linestyle', '--', 'color','w',  'linewidth',3,'HandleVisibility','off')
    line([int{i}{2}(1) int{i}{2}(1)], [int{i}{1}(1) int{i}{1}(2)], 'linestyle', '--', 'color','w',  'linewidth',3,'HandleVisibility','off')
    line([int{i}{2}(2) int{i}{2}(2)], [int{i}{1}(1) int{i}{1}(2)], 'linestyle', '--', 'color','w',  'linewidth',3,'HandleVisibility','off')
end
plot(newY(:,2),newY(:,1), 'wx', 'markersize', 12,'HandleVisibility','off')
colormap(mycmap()); colorbar
xlabel('Second parameter'); ylabel('First parameter')
xlim([intt_(floor(lw/2)+1) intt_(length(intt_)-floor(lw/2))])
ylim([intt_(floor(lw/2)+1) intt_(length(intt_)-floor(lw/2))])
title('(supp)')
    
set(gca, 'fontsize', 18)
set(gca,'DataAspectRatio',[10 10 10])
set(gca,'YDir','normal')

clear mycmap;
mycmap = mycmap();

bounds2 = [-0.2 0.2];

h(9) = subplot(349);
err     = Rmse_grid_ext - Rmse_grid;
imagesc(inter1,inter2,mean(err,3), bounds2)
hold on
for i = 1:length(int)
    line([int{i}{2}(1) int{i}{2}(2)], [int{i}{1}(1) int{i}{1}(1)], 'linestyle', '--', 'color','w',  'linewidth',3,'HandleVisibility','off')
    line([int{i}{2}(1) int{i}{2}(2)], [int{i}{1}(2) int{i}{1}(2)], 'linestyle', '--', 'color','w',  'linewidth',3,'HandleVisibility','off')
    line([int{i}{2}(1) int{i}{2}(1)], [int{i}{1}(1) int{i}{1}(2)], 'linestyle', '--', 'color','w',  'linewidth',3,'HandleVisibility','off')
    line([int{i}{2}(2) int{i}{2}(2)], [int{i}{1}(1) int{i}{1}(2)], 'linestyle', '--', 'color','w',  'linewidth',3,'HandleVisibility','off')
end
plot(newY(:,2),newY(:,1), 'wx', 'markersize', 12,'HandleVisibility','off')
colormap(mycmap); colorbar
xlabel('Second parameter'); ylabel('First parameter')
xlim([intt_(floor(lw/2)+1) intt_(length(intt_)-floor(lw/2))])
ylim([intt_(floor(lw/2)+1) intt_(length(intt_)-floor(lw/2))])
title('Diff')

set(gca, 'fontsize', 18)
set(gca,'DataAspectRatio',[10 10 10])
set(gca,'YDir','normal')


h(10) = subplot(3,4,10);
err     = Rmse_gllim_ext - Rmse_gllim;
imagesc(inter1,inter2,mean(err,3), bounds2)
hold on
for i = 1:length(int)
    line([int{i}{2}(1) int{i}{2}(2)], [int{i}{1}(1) int{i}{1}(1)], 'linestyle', '--', 'color','w',  'linewidth',3,'HandleVisibility','off')
    line([int{i}{2}(1) int{i}{2}(2)], [int{i}{1}(2) int{i}{1}(2)], 'linestyle', '--', 'color','w',  'linewidth',3,'HandleVisibility','off')
    line([int{i}{2}(1) int{i}{2}(1)], [int{i}{1}(1) int{i}{1}(2)], 'linestyle', '--', 'color','w',  'linewidth',3,'HandleVisibility','off')
    line([int{i}{2}(2) int{i}{2}(2)], [int{i}{1}(1) int{i}{1}(2)], 'linestyle', '--', 'color','w',  'linewidth',3,'HandleVisibility','off')
end
plot(newY(:,2),newY(:,1), 'wx', 'markersize', 12,'HandleVisibility','off')
colormap(mycmap); colorbar
xlabel('Second parameter'); ylabel('First parameter')
xlim([intt_(floor(lw/2)+1) intt_(length(intt_)-floor(lw/2))])
ylim([intt_(floor(lw/2)+1) intt_(length(intt_)-floor(lw/2))])
title('Diff')

set(gca, 'fontsize', 18)
set(gca,'DataAspectRatio',[10 10 10])
set(gca,'YDir','normal')


h(11) = subplot(3,4,11);
err     = Rmse_nn_ext - Rmse_nn;
imagesc(inter1,inter2,mean(err,3), bounds2)
hold on
for i = 1:length(int)
    line([int{i}{2}(1) int{i}{2}(2)], [int{i}{1}(1) int{i}{1}(1)], 'linestyle', '--', 'color','w',  'linewidth',3,'HandleVisibility','off')
    line([int{i}{2}(1) int{i}{2}(2)], [int{i}{1}(2) int{i}{1}(2)], 'linestyle', '--', 'color','w',  'linewidth',3,'HandleVisibility','off')
    line([int{i}{2}(1) int{i}{2}(1)], [int{i}{1}(1) int{i}{1}(2)], 'linestyle', '--', 'color','w',  'linewidth',3,'HandleVisibility','off')
    line([int{i}{2}(2) int{i}{2}(2)], [int{i}{1}(1) int{i}{1}(2)], 'linestyle', '--', 'color','w',  'linewidth',3,'HandleVisibility','off')
end
plot(newY(:,2),newY(:,1), 'wx', 'markersize', 12,'HandleVisibility','off')
colormap(mycmap); colorbar
xlabel('Second parameter'); ylabel('First parameter')
xlim([intt_(floor(lw/2)+1) intt_(length(intt_)-floor(lw/2))])
ylim([intt_(floor(lw/2)+1) intt_(length(intt_)-floor(lw/2))])
title('Diff')

set(gca, 'fontsize', 18)
set(gca,'DataAspectRatio',[10 10 10])
set(gca,'YDir','normal')


 
%%

if backup == 1
    savefig(fig, ['figures/' 'BoundaryBehavior'])
end


%%

% figure;
% subplot(121);
% plot(reshape(mah,1,[]), reshape(mean(err,3),1,[]), '.')
% subplot(122);
% plot(reshape(mah,1,[]), reshape(mean(CI(:,:,1),3),1,[]), '.')
% 
% if length(numel(mask)) == 3
%     mask = mask(:,:,1) + mask(:,:,2);
%     mask_out = mask_out(:,:,1) & mask_out(:,:,2);
% end
% 
% msk = mask_out;
% 
% err_patches = (mean(Rmse_grid,3) .* msk);
% err_patches(err_patches == 0) = nan;
% err_patches_ext = (mean(Rmse_grid_ext,3) .* msk);
% err_patches_ext(err_patches_ext == 0) = nan;
% 
% disp('Relat error (%) with DM:')
% disp(100* (nanmean(err_patches(:)) - nanmean(err_patches_ext(:))) ./nanmean(err_patches(:)))
% 
% err_patches = (mean(Rmse_nn,3) .* msk);
% err_patches(err_patches == 0) = nan;
% err_patches_ext = (mean(Rmse_nn_ext,3) .* msk);
% err_patches_ext(err_patches_ext == 0) = nan;
% 
% disp('Relat error (%) with DB-DL:')
% disp(100* (nanmean(err_patches(:)) - nanmean(err_patches_ext(:))) ./nanmean(err_patches(:)))
% 
% err_patches = (mean(Rmse_gllim,3) .* msk);
% err_patches(err_patches == 0) = nan;
% err_patches_ext = (mean(Rmse_gllim_ext,3) .* msk);
% err_patches_ext(err_patches_ext == 0) = nan;
% 
% disp('Relat error (%) with DB-SL:')
% disp(100* (nanmean(err_patches(:)) - nanmean(err_patches_ext(:))) ./nanmean(err_patches(:)))


