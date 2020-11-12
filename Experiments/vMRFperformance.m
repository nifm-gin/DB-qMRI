
%% Description
%
% Blabla
%
% Fabien Boux - 11/2020

Init
disp(['Running experiment ' mfilename '.m'])


%% Parameters

verbose = 0;
backup = 1;

% Function definition
nb_param = 3;
nb_train_signals = 150500;
nb_test_signals  = 100000;

snr_test  	= 100;
snr_train   = 60;

dicos   = {{'inputs/new/qMC.mat', ' '},... %qMC
           {'inputs/new/Regular.mat', ' '},... %grid
           {'inputs/new/qMC.mat', ' '},... %test signals
          };
use_recomputed_parameter = 0;
     
int{1}  = [-1.8 28] *1e-2;
int{2}  = [-1.8 48] *1e-6;
N       = 150;
lw      = 16;

vsi_bounds = {[2   7]*1e-6, ... %for bvf estimation
              [7  15]*1e-6, ...
              [15 30]*1e-6};
          
bvf_bounds = {[.5  4]*1e-2, ... %for vsi estimation
              [4  10]*1e-2, ...
              [10 20]*1e-2};


%%

Rmse_gllim = nan(N,length(vsi_bounds),length(dicos{1}),nb_param);
Rmse_grid   = Rmse_gllim; Rmse_anlt = Rmse_gllim; Rmse_nn = Rmse_gllim; 
Mae_gllim  = Rmse_gllim; Mae_grid = Rmse_gllim; Mae_anlt = Rmse_gllim; Mae_nn = Rmse_gllim;
CI         = Rmse_gllim;

for s = 1:1
    
    % Load MGEFIDSE signals
    
    % Load echotimes
    Te     	= load('inputs/echoetime.mat');
    Te    	= Te.t;
    Tse     = 0.060; %sec
    
    % Quasi-random
    load(dicos{1}{s})
    if use_recomputed_parameter == 1
        Dico.Parameters.Par(:,1) = Dico.Parameters.Par(:,7);
        Dico.Parameters.Par(:,2) = Dico.Parameters.Par(:,8);
    end
    Xqmc    = Dico.MRSignals(:,end/2+1:end) ./ Dico.MRSignals(:,1:end/2);
    parfor i = 1:size(Xqmc,1)
        Tmp(i,:)    = interp1(Dico.Tacq(1:size(Xqmc,2)), Xqmc(i,:), Te);
    end
    Xqmc    = Tmp; Tmp  = [];
    Yqmc    = Dico.Parameters.Par;

    [rows, ~] = find(~isnan(Xqmc));
    r       = unique(rows);
    Xqmc    = Xqmc(r,:);
    Yqmc    = Yqmc(r,1:nb_param);
    
    r       = unique(find(any(isnan(Yqmc),2) == 0));
    Xqmc    = Xqmc(r,:);
    Yqmc    = Yqmc(r,1:nb_param);

    Yqmc  = Yqmc(:, [2 3 1]);
    
    % Regular
    load(dicos{2}{s})
    if use_recomputed_parameter == 1
        Dico.Parameters.Par(:,1) = Dico.Parameters.Par(:,7);
        Dico.Parameters.Par(:,2) = Dico.Parameters.Par(:,8);
    end
    Xgrid 	= Dico.MRSignals(:,end/2+1:end) ./ Dico.MRSignals(:,1:end/2);
    parfor i = 1:size(Xgrid,1)
        Tmp(i,:)    = interp1(Dico.Tacq(1:size(Xgrid,2)), Xgrid(i,:), Te);
    end
    Xgrid   = Tmp; Tmp  = [];
    Ygrid   = Dico.Parameters.Par;

    [rows, ~] = find(~isnan(Xgrid));
    r       = unique(rows);
    Xgrid   = Xgrid(r,:);
    Ygrid   = Ygrid(r,1:nb_param);
    
    r       = unique(find(any(isnan(Ygrid),2) == 0));
    Xgrid   = Xgrid(r,:);
    Ygrid   = Ygrid(r,1:nb_param);
    
    Ygrid  = Ygrid(:, [2 3 1]);
    
    % Picking test signals
    load(dicos{3}{s})
    if use_recomputed_parameter == 1
        Dico.Parameters.Par(:,1) = Dico.Parameters.Par(:,7);
        Dico.Parameters.Par(:,2) = Dico.Parameters.Par(:,8);
    end
    Xtest_  = Dico.MRSignals;
    parfor i = 1:size(Xtest_,1)
        Tmp1(i,:)    = interp1(Dico.Tacq(1:size(Xtest_,2)/2), Xtest_(i,1:end/2), Te);
        Tmp2(i,:)    = interp1(Dico.Tacq(1:size(Xtest_,2)/2), Xtest_(i,end/2+1:end), Te);
    end
    Xtest_  = [Tmp1 Tmp2]; Tmp1 = []; Tmp2 = [];
    Ytest_  = Dico.Parameters.Par;
    
    [rows, ~] = find(~isnan(Xtest_));
    r       = unique(rows);    
    Xtest_  = Xtest_(r,:);
    Ytest_  = Ytest_(r,1:nb_param);
    
    Xtest_ = Xtest_(1:nb_test_signals,:);
    Ytest_ = Ytest_(1:nb_test_signals,[2 3 1]);
        
    
    %% DBL method training
    
    Model_.K    = 50;
    Model_.Lw   = 0;
    
    r           = randperm(size(Xqmc,1));
    Dico        = FormatDico(Xqmc(r(1:nb_train_signals),:), Yqmc(r(1:nb_train_signals),:));
    [~,Model]   = AnalyzeMRImages([], Dico, 'DB-SL', Model_);
    
    [~,NeuralNet] = AnalyzeMRImages([], Dico, 'DB-DL');
    
    %prepare dico for matching method
    Dico = FormatDico(Xgrid(1:nb_train_signals,:), Ygrid(1:nb_train_signals,:));
    

    %% Estimation of BVf

    for c = 1:length(vsi_bounds) 

        int_	= int{1}(1): (int{1}(2) - int{1}(1))/(N-1) :int{1}(2);

        parfor m = 1+floor(lw/2):length(int_)-floor(lw/2)

            disp([1 c m])

            subint  = int_([m-floor(lw/2)  m+floor(lw/2)-1]);
            r       = (Ytest_(:,1) >= subint(1))     & (Ytest_(:,1) <= subint(2))... % BVf condition
                      & (Ytest_(:,2) >= vsi_bounds{c}(1)) & (Ytest_(:,2) <= vsi_bounds{c}(2)); % VSI condition

            dens(m,c,s,1) = sum(r);

            if sum(r) > 1
                Ytest   = Ytest_(r,:);
                Xtest   = AddNoise(Xtest_(r,:), snr_test);
                
                %CEF
                Y       = EstimateParametersFromAnalytic(Te,Tse, Xtest(:,1:end/2), Xtest(:,end/2+1:end), []);
                Y(:,1)  = Y(:,1) * 1e3;
                Y( Y(:,1) < 0 | Y(:,2) < 0 | Y(:,1) > int{1}(2) | Y(:,2) > int{2}(2), :) = nan;
                [Rmse_anlt(m,c,s,1),~, Mae_anlt(m,c,s,1)]  = EvaluateEstimation(Ytest(:,1), Y(:,1));
                
                %dictionary-based methods (first ratio is computed)
                Xtest   = Xtest(:,end/2+1:end)./Xtest(:,1:end/2);
                
                %DBM
                Estim   = AnalyzeMRImages(Xtest, Dico, 'DBM');
                Y       = Estim.GridSearch.Y(:,1:nb_param);
                [Rmse_grid(m,c,s,1),~, Mae_grid(m,c,s,1)]  = EvaluateEstimation(Ytest(:,1), Y(:,1));
                
                %DB-SL
                Estim   = AnalyzeMRImages(Xtest, [], 'DB-SL', Model, [],[], snr_test); 
                Y       = Estim.Regression.Y(:,1:nb_param);
                [Rmse_gllim(m,c,s,1),~, Mae_gllim(m,c,s,1)] = EvaluateEstimation(Ytest(:,1), Y(:,1));
                CI(m,c,s,1) 	= nanmean(Estim.Regression.Cov(:,1,1)).^.5;            

                %DB-DL
                Estim   = AnalyzeMRImages(Xtest, [], 'DB-DL', NeuralNet); 
                Y       = Estim.Regression.Y(:,1:nb_param);
                [Rmse_nn(m,c,s,1),~, Mae_nn(m,c,s,1)] = EvaluateEstimation(Ytest(:,1), Y(:,1));
            end

        end %BVf sliding window
    end %VSI cases


    %% Estimation of VSI
    for c = 1:length(bvf_bounds) 

        int_	= int{2}(1): (int{2}(2) - int{2}(1))/(N-1) :int{2}(2);

        parfor m = 1+floor(lw/2):length(int_)-floor(lw/2)

            disp([2 c m])

            subint  = int_([m-floor(lw/2)  m+floor(lw/2)-1]);
            r       = (Ytest_(:,2) >= subint(1))     & (Ytest_(:,2) <= subint(2))... % BVf condition
                      & (Ytest_(:,1) >= bvf_bounds{c}(1)) & (Ytest_(:,1) <= bvf_bounds{c}(2)); % VSI condition

            dens(m,c,s,2) = sum(r);

            if sum(r) > 1
                Ytest   = Ytest_(r,:);
                Xtest   = AddNoise(Xtest_(r,:), snr_test);

                %CEF
                Y       = EstimateParametersFromAnalytic(Te,Tse, Xtest(:,1:end/2), Xtest(:,end/2+1:end), []);
                Y(:,1)  = Y(:,1) * 1e3;
                Y( Y(:,1) < 0 | Y(:,2) < 0 | Y(:,1) > int{1}(2) | Y(:,2) > int{2}(2), :) = nan;
                [Rmse_anlt(m,c,s,2),~, Mae_anlt(m,c,s,2)]  = EvaluateEstimation(Ytest(:,2), Y(:,2));
                
                %dictionary-based methods (first ratio is computed)
                Xtest   = Xtest(:,end/2+1:end)./Xtest(:,1:end/2);
                
                %DBM
                Estim   = AnalyzeMRImages(Xtest, Dico, 'DBM');
                Y       = Estim.GridSearch.Y(:,1:nb_param);
                [Rmse_grid(m,c,s,2),~, Mae_grid(m,c,s,2)]  = EvaluateEstimation(Ytest(:,2), Y(:,2));
                
                %DB-SL
                Estim   = AnalyzeMRImages(Xtest, [], 'DB-SL', Model, [],[], snr_test); 
                Y       = Estim.Regression.Y(:,1:nb_param);
                [Rmse_gllim(m,c,s,2),~, Mae_gllim(m,c,s,2)] = EvaluateEstimation(Ytest(:,2), Y(:,2));
                CI(m,c,s,2) 	= nanmean(Estim.Regression.Cov(:,1,2)).^.5;            

                %DB-DL
                Estim   = AnalyzeMRImages(Xtest, [], 'DB-DL', NeuralNet); 
                Y       = Estim.Regression.Y(:,1:nb_param);
                [Rmse_nn(m,c,s,2),~, Mae_nn(m,c,s,2)] = EvaluateEstimation(Ytest(:,2), Y(:,2));
            end

        end %VSI sliding window
    end %BVf cases
end %StO2 cases


%% Saving 

if backup == 1
    clear tmp* Dico* Estim X* Y* Rmse_ Mae_
    save(['temp/' 'vMRFperformance'])
end


%% Plot

tit{1}  = {'Small vessels', 'Medium vessels', 'Large vessels'};
tit{2}  = {'Low BVf', 'Average BVf', 'High BVf'};
xl      = {'BVf (%)'};
yl      = {'RMSE on BVf (%)', 'RMSE on VSI (µm)'};

% Remove unexpected values
vsi_limits = {[2  50]*1e-6,...
              [4  50]*1e-6,...
              [4  50]*1e-6};
bvf_limits = {[.5 25]*1e-2,...
              [.5 25]*1e-2,...
              [.5 25]*1e-2};

colors  = [          0    0.4470    0.7410
                0.8500    0.3250    0.0980
                0.9290    0.6940    0.1250
                0.4940    0.1840    0.5560
                0.4660    0.6740    0.1880
                0.3010    0.7450    0.9330
                0.6350    0.0780    0.1840];
colors  = [colors; colors];


fig = figure;
for c = 1:length(vsi_limits)
        
    int_	= int{1}(1): (int{1}(2) - int{1}(1))/(N-1) :int{1}(2);

    h(c) = subplot(2,length(vsi_limits), c);
    hold on
    plot(1e2 * int_, 1e2 * Rmse_anlt(:,c,s,1),  '-', 'linewidth',2, 'color', colors(1,:))
    plot(1e2 * int_, 1e2 * Rmse_grid(:,c,s,1),  '-', 'linewidth',2, 'color', colors(2,:))
    plot(1e2 * int_, 1e2 * Rmse_nn(:,c,s,1),    '-', 'linewidth',2, 'color', colors(4,:))
    plot(1e2 * int_, 1e2 * Rmse_gllim(:,c,s,1), '-', 'linewidth',2, 'color', colors(3,:))
    
    plot(1e2 * int_, 1e2 * CI(:,c,s,1), '--', 'linewidth',2, 'color', colors(3,:))    
    
    xlabel('BVf (%)')

    if s == 1
        title(tit{1}{c})
    end
    if s == length(dicos{1})
        xlabel(xl{1})
    end
    if mod(c,(length(dicos{1})-1)*(length(vsi_limits))) == 1
        ylabel(yl{1})
    end        

    switch c
        case 1
            xlim([1 21])
        case 2
            xlim([1 21])
        case 3
            xlim([1 21])
    end

    if c == length(vsi_limits) && s == 1
        legend('CEF','DBM','DB-DL','DB-SL')
    end

end %VSI cases

for c = 1:length(bvf_limits)
        
    int_	= int{2}(1): (int{2}(2) - int{2}(1))/(N-1) :int{2}(2);

    h(length(bvf_limits)+c) = subplot(2,length(bvf_limits), length(bvf_limits)+c);
    hold on
    plot(1e6 * int_, 1e6 * Rmse_anlt(:,c,s,2),  '-', 'linewidth',2, 'color', colors(1,:))
    plot(1e6 * int_, 1e6 * Rmse_grid(:,c,s,2),  '-', 'linewidth',2, 'color', colors(2,:))
    plot(1e6 * int_, 1e6 * Rmse_nn(:,c,s,2),    '-', 'linewidth',2, 'color', colors(4,:))
    plot(1e6 * int_, 1e6 * Rmse_gllim(:,c,s,2), '-', 'linewidth',2, 'color', colors(3,:))
    
    plot(1e6 * int_, 1e6 * CI(:,c,s,2), '--', 'linewidth',2, 'color', colors(3,:))   
    
    xlabel('VSI (µm)')
    
    if s == 1
        title(tit{2}{c})
    end
    if s == length(dicos{1})
        xlabel(xl{1})
    end
    if mod(c,(length(dicos{1})-1)*(length(vsi_limits))) == 1
        ylabel(yl{2})
    end        

    switch c 
        case 1
            xlim([1 40])
        case 2
            xlim([1 40])
        case 3
            xlim([1 40])
    end
end %BVf cases

linkaxes(h(1:3),'xy')
linkaxes(h(4:6),'xy')
set(h, 'fontsize',15)


%% Exporting figures

if backup == 1
    savefig(fig, ['figures/' 'vMRFperformance'])
end
