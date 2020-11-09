
%% Description
%
% We repeat the experiment of BoundaryBehavior.m script (see this file for 
% complete description) using vascular MRF signals (ratio MGEFIDSE pre/post
% USPIO injection).
%
% Fabien Boux - 11/2020

Init
disp(['Running experiment ' mfilename '.m'])


%% Parameters

backup  = 1;
verbose = 1;
display = 0; %if 1 plot extra figures, else only the main figure

%experiment settings
dbdl_computation = 1;
use_extra_entries = 1;
N       = 150;
lw      = 8;
snr     = 60;
nb_param = 3;

%regression settings
Model_.K    = 60;
Model_.Lw   = 0;


%% Load MGEFIDSE signals

%dictionary filenames
dicos   = {'inputs/new/qMC.mat',... %qMC
           'inputs/new/Regular.mat',... %grid
           'inputs/new/qMC.mat',... %test signals
          };
      
% Load echotimes
echtime = load('inputs/echoetime.mat');
echtime = echtime.t;

% Quasi-random
load(dicos{1})
Xqmc    = Dico.MRSignals(:,end/2+1:end) ./ Dico.MRSignals(:,1:end/2);
% parfor i = 1:size(Xqmc,1)
%     Tmp(i,:)    = interp1(Dico.Tacq(1:size(Xqmc,2)), Xqmc(i,:), echtime);
% end
% Xqmc    = Tmp; Tmp  = [];
Yqmc    = Dico.Parameters.Par;

Xqmc(Dico.Parameters.Par(:,4)== 1e-6,:) = [];
Yqmc(Dico.Parameters.Par(:,4)== 1e-6,:) = [];

Yqmc  = Yqmc(:, [2 3 1]);

[rows, ~] = find(~isnan(Xqmc));
r       = unique(rows);
Xqmc    = Xqmc(r,:);
Yqmc    = Yqmc(r,:);

r       = unique(find(any(isnan(Yqmc),2) == 0));
Xqmc    = Xqmc(r,:);
Yqmc    = Yqmc(r,:);


% Grid
load(dicos{2})
Xgrid 	= Dico.MRSignals(:,end/2+1:end) ./ Dico.MRSignals(:,1:end/2);
% parfor i = 1:size(Xgrid,1)
%     Tmp(i,:)    = interp1(Dico.Tacq(1:size(Xgrid,2)), Xgrid(i,:), echtime);
% end
% Xgrid   = Tmp; Tmp  = [];
Ygrid   = Dico.Parameters.Par;

Xgrid(Dico.Parameters.Par(:,4)== 1e-6,:) = [];
Ygrid(Dico.Parameters.Par(:,4)== 1e-6,:) = [];

Ygrid  = Ygrid(:, [2 3 1]);

[rows, ~] = find(~isnan(Xgrid));
r       = unique(rows);
Xgrid   = Xgrid(r,:);
Ygrid   = Ygrid(r,:);

r       = unique(find(any(isnan(Ygrid),2) == 0));
Xgrid   = Xgrid(r,:);
Ygrid   = Ygrid(r,:);


% Random - test signals
load(dicos{3})
Xtest_ 	= Dico.MRSignals(:,end/2+1:end) ./ Dico.MRSignals(:,1:end/2);
% parfor i = 1:size(Xtest_,1)
%     Tmp(i,:)    = interp1(Dico.Tacq(1:size(Xtest_,2)), Xtest_(i,:), echtime);
% end
% Xtest_  = Tmp; Tmp  = [];
Ytest_  = Dico.Parameters.Par;

Xtest_(Dico.Parameters.Par(:,4)== 1e-6,:) = [];
Ytest_(Dico.Parameters.Par(:,4)== 1e-6,:) = [];

[rows, ~] = find(~isnan(Xtest_));
r       = unique(rows);    
Xtest_  = Xtest_(r,:);
Ytest_  = Ytest_(r,:);

Ytest_ = Ytest_(:,[2 3 1]);

Xtest_(Ytest_(:,2) > 40 *1e-6,:) = [];
Ytest_(Ytest_(:,2) > 40 *1e-6,:) = [];


%% Reduce datasets to keep only signals with StO2 around 0.7

bds_sto2 = [0.65 0.75];
bds_sto2 = [0.6 0.8];

v   = (Ygrid(:,3) <= bds_sto2(1) | Ygrid(:,3) >= bds_sto2(2));
Xgrid(v,:)  = [];
Ygrid(v,:)  = [];

v   = (Yqmc(:,3) <= bds_sto2(1) | Yqmc(:,3) >= bds_sto2(2));
Xqmc(v,:)   = [];
Yqmc(v,:)   = [];

v   = (Ytest_(:,3) <= bds_sto2(1) | Ytest_(:,3) >= bds_sto2(2));
Xtest_(v,:) = [];
Ytest_(v,:) = [];


%% Reduce dictionaries to 2 patches & Add extra entries if required

int_bvf = {[0.03 0.07], [0.15 0.25]};
int_vsi = {[4e-6 12e-6], [12e-6 20e-6]};

%extra entries
if use_extra_entries == 1
            %some entries in space around patches
    ext = [ 20e-2    5e-6     ;
            12e-2    2.7e-6   ;
            10e-2    9e-6     ;
             4e-2   20e-6     ;
             7e-2   21e-6     ;
            10e-2   22e-6     ;
            13e-2   23e-6     ;
            16e-2   24e-6     ;
            19e-2   25e-6     ;
            22e-2   26e-6     ;
            ... %then, boundary entries
            max(Yqmc(:,1)) 5e-6;
            max(Yqmc(:,1)) max(Yqmc(:,2));
            15e-2 max(Yqmc(:,2));
            ];
else
    ext = [];
end

if  ~isempty(ext) == 1
    for i = 1:size(ext,1)
        Xtmp = [Xtest_; Xqmc];
        Ytmp = [Ytest_(:,1:2); Yqmc(:,1:2)];
        Yttmp = [Ytest_; Yqmc];
        mi = min(Ytmp); ma = max(Ytmp - mi);
        Ytmp = (Ytmp - mi) ./ ma;
        lg = find(sum((Ytmp - ((ext(i,:)-mi) ./ma)).^2,2) == min(sum((Ytmp - ((ext(i,:)-mi) ./ma)).^2,2)));

        Xgrid_ext(i,:) = Xtmp(lg(ceil(end/2)),:);
        Ygrid_ext(i,:) = Yttmp(lg(ceil(end/2)),:);
        Xqmc_ext(i,:)  = Xtmp(lg(ceil(end/2)),:);
        Yqmc_ext(i,:)  = Yttmp(lg(ceil(end/2)),:);
    end
else
    Xgrid_ext = [];
    Ygrid_ext = [];
    Xqmc_ext  = [];
    Yqmc_ext  = [];
end

%patches
for i = 1:numel(int_bvf)
    v_grid(:,i)   = (( int_bvf{i}(1) <= Ygrid(:,1) ) &  ( Ygrid(:,1) <= int_bvf{i}(2) )) & ...
                    (( int_vsi{i}(1) <= Ygrid(:,2) ) &  ( Ygrid(:,2) <= int_vsi{i}(2) ));
    v_qmc(:,i)    = (( int_bvf{i}(1) <= Yqmc(:,1) )  &  ( Yqmc(:,1) <= int_bvf{i}(2) )) & ...
                    (( int_vsi{i}(1) <= Yqmc(:,2) )  &  ( Yqmc(:,2) <= int_vsi{i}(2) ));
end

Xgrid(~any(v_grid'),:) = [];
Ygrid(~any(v_grid'),:) = [];
Xqmc(~any(v_qmc'),:)  = [];
Yqmc(~any(v_qmc'),:)  = [];

%concat patches and extra entries
Xgrid = [Xgrid; Xgrid_ext];
Ygrid = [Ygrid; Ygrid_ext];
Xqmc  = [Xqmc;  Xqmc_ext];
Yqmc  = [Yqmc;  Yqmc_ext];

if display == 1
    figure
    hold on
    plot(Ygrid(:,1), Ygrid(:,2), '.')
    for i = 1:numel(int_bvf)
        plot([int_bvf{i}(1) int_bvf{i}(2)], [int_vsi{i}(1) int_vsi{i}(1)], 'k--', 'linewidth',2)
        plot([int_bvf{i}(1) int_bvf{i}(2)], [int_vsi{i}(2) int_vsi{i}(2)], 'k--', 'linewidth',2)
        plot([int_bvf{i}(1) int_bvf{i}(1)], [int_vsi{i}(1) int_vsi{i}(2)], 'k--', 'linewidth',2)
        plot([int_bvf{i}(2) int_bvf{i}(2)], [int_vsi{i}(1) int_vsi{i}(2)], 'k--', 'linewidth',2)
    end
    if ~isempty(ext) == 1
        plot(Yqmc_ext(:,1),Yqmc_ext(:,2), 'x')
    end
end
clear *tmp* 


%% DBL methods training

Dico    = FormatDico(Xqmc, Yqmc(:,1:nb_param));

%DB-SL
[~,Model] = AnalyzeMRImages([], Dico, 'DB-SL', Model_);

%DB-DL
if dbdl_computation == 1
    [~,NeuralNet] = AnalyzeMRImages([], Dico, 'DB-DL');
end


%% Estimations

for i = 1:2
    intt_{i} = min(Ytest_(:,i)) : (max(Ytest_(:,i)) - min(Ytest_(:,i))) / N : max(Ytest_(:,i));
end

%init
dens        = zeros(length(intt_{1})-floor(lw/2), length(intt_{2})-floor(lw/2));
inter1      = zeros(1,size(dens,1));     
inter2      = zeros(1,size(dens,2));
Rmse_grid   = nan(size(dens,1), size(dens,2), nb_param);
Rmse_gllim  = Rmse_grid;
Rmse_nn     = Rmse_grid;
Mae_grid    = Rmse_grid; 
Mae_gllim   = Rmse_grid;
Mae_nn      = Rmse_grid;
CI          = Rmse_grid;

%run
for s1 = floor(lw/2)+1:length(intt_{1})-floor(lw/2)
    
    if verbose == 1, disp([num2str(s1-floor(lw/2)) '/' num2str(length(floor(lw/2)+1:length(intt_{1})-floor(lw/2)))]); end

    inter1(s1) = intt_{1}(s1);
    
    subint1 = intt_{1}([s1-floor(lw/2) s1+floor(lw/2)]);
    
    v1      = (subint1(1) <= Ytest_(:,1)) & (Ytest_(:,1) < subint1(2));
    Xtest__ = Xtest_(v1,:);
    Ytest__ = Ytest_(v1,:);
    
    for s2 = floor(lw/2)+1:length(intt_{2})-floor(lw/2)

%         % Test data
        subint2 = intt_{2}([s2-floor(lw/2) s2+floor(lw/2)]);
        inter2(s2) = intt_{2}(s2);

        v2  	= (subint2(1) <= Ytest__(:,2)) & (Ytest__(:,2) < subint2(2));
        
        Xtest   = AddNoise(Xtest__(v2,:), snr);
        Ytest   = Ytest__(v2,:);
        dens(s1,s2) = sum(v2(:));
        
        for i = 1:2
            mask(s1,s2,i)   = (( int_bvf{i}(1) <= subint1(1) ) &  ( subint1(2) <= int_bvf{i}(2) )) & ...
                            (( int_vsi{i}(1) <= subint2(1) ) &  ( subint2(2) <= int_vsi{i}(2) ));
            mask_neighboor(s1,s2,i) = (( (int_bvf{i}(1) -0.02) <= subint1(1) ) &  ( subint1(2) <= (int_bvf{i}(2)+0.02) )) & ...
                            (( (int_vsi{i}(1) -2e-6)  <= subint2(1) ) &  ( subint2(2) <= (int_vsi{i}(2) +2e-6)));
            mask_dist(s1,s2,i) = (( (int_bvf{i}(1) -0.02) <= subint1(1) ) &  ( subint1(2) <= (int_bvf{i}(2)+0.02) )) & ...
                            (( (int_vsi{i}(1) -2e-6)  <= subint2(1) ) &  ( subint2(2) <= (int_vsi{i}(2) +2e-6)));
        end
        
        if dens(s1,s2) > 1
            % Compute estimates
            %DBM
            Y       = EstimateParametersFromGrid(Xtest,Xgrid,Ygrid);
            [Rmse_grid(s1,s2,:),~, Mae_grid(s1,s2,:)] = EvaluateEstimation(Ytest, Y);
            
            %DB-sL
            Estim   = AnalyzeMRImages(Xtest, [], 'DB-SL', Model, [],[], snr);
            [Rmse_gllim(s1,s2,:),~, Mae_gllim(s1,s2,:)] = EvaluateEstimation(Ytest(:,1:nb_param), squeeze(Estim.Regression.Y));
            CI(s1,s2,:) = nanmean(squeeze(Estim.Regression.Cov.^.5));

            %DB-DL
            if dbdl_computation == 1
                Estim   = AnalyzeMRImages(Xtest, [], 'DB-DL', NeuralNet);
                [Rmse_nn(s1,s2,:),~, Mae_nn(s1,s2,:)] = EvaluateEstimation(Ytest(:,1:nb_param), squeeze(Estim.Regression.Y));
            end
        end
    end
end
    

%% Saving 

if backup == 1
    clear tmp* X* Dic
    if isempty(ext)
        save(['temp/' 'BoundaryBehaviour_vMRF_1'])
    else
        save(['temp/' 'BoundaryBehaviour_vMRF_2'])
    end
end


%% Display

fig = figure;
mycmap = mycmap_extended;

for p = 1:2
    
    switch p
        case 1
            bounds = [0 15] *1e-2;
        case 2
            bounds = [0 20] * 1e-6;
    end
    
    %DBM
    h(4*(p-1)+1) = subplot(2,4,4*(p-1)+1);
    err     = Rmse_grid(:,:,p);

    imagesc(inter1, inter2, err', bounds)
    hold on
    for i = 1:length(int_bvf)
        line([int_bvf{i}(1) int_bvf{i}(2)], [int_vsi{i}(1) int_vsi{i}(1)], 'linestyle', '--', 'color','w',  'linewidth',3,'HandleVisibility','off')
        line([int_bvf{i}(1) int_bvf{i}(2)], [int_vsi{i}(2) int_vsi{i}(2)], 'linestyle', '--', 'color','w',  'linewidth',3,'HandleVisibility','off')
        line([int_bvf{i}(1) int_bvf{i}(1)], [int_vsi{i}(1) int_vsi{i}(2)], 'linestyle', '--', 'color','w',  'linewidth',3,'HandleVisibility','off')
        line([int_bvf{i}(2) int_bvf{i}(2)], [int_vsi{i}(1) int_vsi{i}(2)], 'linestyle', '--', 'color','w',  'linewidth',3,'HandleVisibility','off')
    end 
    if ~isempty(Ygrid_ext)
        plot(Ygrid_ext(:,1),Ygrid_ext(:,2), 'wx')
    end
    colormap(mycmap()); colorbar
    xlabel('BVf (%)'); ylabel('VSI (um)')
    title('(a) DBM')

    
    %DB-SL
    h(4*(p-1)+2) = subplot(2,4,4*(p-1)+2);
    err     = Rmse_gllim(:,:,p);
    
    imagesc(inter1, inter2, err', bounds)
    hold on
    for i = 1:length(int_bvf)
        line([int_bvf{i}(1) int_bvf{i}(2)], [int_vsi{i}(1) int_vsi{i}(1)], 'linestyle', '--', 'color','w',  'linewidth',3,'HandleVisibility','off')
        line([int_bvf{i}(1) int_bvf{i}(2)], [int_vsi{i}(2) int_vsi{i}(2)], 'linestyle', '--', 'color','w',  'linewidth',3,'HandleVisibility','off')
        line([int_bvf{i}(1) int_bvf{i}(1)], [int_vsi{i}(1) int_vsi{i}(2)], 'linestyle', '--', 'color','w',  'linewidth',3,'HandleVisibility','off')
        line([int_bvf{i}(2) int_bvf{i}(2)], [int_vsi{i}(1) int_vsi{i}(2)], 'linestyle', '--', 'color','w',  'linewidth',3,'HandleVisibility','off')
    end 
    if ~isempty(Yqmc_ext)
        plot(Yqmc_ext(:,1),Yqmc_ext(:,2), 'wx')
    end
    colormap(mycmap()); colorbar
    xlabel('BVf (%)'); ylabel('VSI (um)')
    title('(b) DB-SL')

    
    %DB-DL
	if dbdl_computation == 1
        h(4*(p-1)+3) = subplot(2,4,4*(p-1)+3);
        err     = Rmse_nn(:,:,p);

        imagesc(inter1, inter2, err', bounds)
        hold on
        for i = 1:length(int_bvf)
            line([int_bvf{i}(1) int_bvf{i}(2)], [int_vsi{i}(1) int_vsi{i}(1)], 'linestyle', '--', 'color','w',  'linewidth',3,'HandleVisibility','off')
            line([int_bvf{i}(1) int_bvf{i}(2)], [int_vsi{i}(2) int_vsi{i}(2)], 'linestyle', '--', 'color','w',  'linewidth',3,'HandleVisibility','off')
            line([int_bvf{i}(1) int_bvf{i}(1)], [int_vsi{i}(1) int_vsi{i}(2)], 'linestyle', '--', 'color','w',  'linewidth',3,'HandleVisibility','off')
            line([int_bvf{i}(2) int_bvf{i}(2)], [int_vsi{i}(1) int_vsi{i}(2)], 'linestyle', '--', 'color','w',  'linewidth',3,'HandleVisibility','off')
        end 
        if ~isempty(Yqmc_ext)
            plot(Yqmc_ext(:,1),Yqmc_ext(:,2), 'wx')
        end
        colormap(mycmap()); colorbar
        xlabel('BVf (%)'); ylabel('VSI (um)')
        title('(c) DB-DL')
    end

    
    h(4*(p-1)+4) = subplot(2,4,4*(p-1)+4);
%     err     = CI(:,:,p);
%     
%     imagesc(inter1, inter2, err', bounds)
%     hold on
%     for i = 1:length(int_bvf)
%         line([int_bvf{i}(1) int_bvf{i}(2)], [int_vsi{i}(1) int_vsi{i}(1)], 'linestyle', '--', 'color','w',  'linewidth',3,'HandleVisibility','off')
%         line([int_bvf{i}(1) int_bvf{i}(2)], [int_vsi{i}(2) int_vsi{i}(2)], 'linestyle', '--', 'color','w',  'linewidth',3,'HandleVisibility','off')
%         line([int_bvf{i}(1) int_bvf{i}(1)], [int_vsi{i}(1) int_vsi{i}(2)], 'linestyle', '--', 'color','w',  'linewidth',3,'HandleVisibility','off')
%         line([int_bvf{i}(2) int_bvf{i}(2)], [int_vsi{i}(1) int_vsi{i}(2)], 'linestyle', '--', 'color','w',  'linewidth',3,'HandleVisibility','off')
%     end 
%     if ~isempty(Yqmc_ext)
%         plot(Yqmc_ext(:,1),Yqmc_ext(:,2), 'wx')
%     end
%     colormap(mycmap()); colorbar
%     xlabel('BVf (%)'); ylabel('VSI (um)')
%     title('(d)')    
    

end
% set(h, 'fontsize', 18)
% set(h,'DataAspectRatio',[10 10 10])
set(h,'YDir','normal')

 
%%

if backup == 1
    if isempty(ext)
        savefig(fig, ['figures/' 'BoundaryBehaviour_vMRF_1'])
    else
        savefig(fig, ['figures/' 'BoundaryBehaviour_vMRF_2'])
    end
end


%% Extra ideas

if numel(size(mask)) == 3
    mask = mask(:,:,1) + mask(:,:,2);
    mask_neighboor = mask_neighboor(:,:,1) + mask_neighboor(:,:,2);
    mask_neighboor = mask_neighboor - mask;
end

[X,Y]   = meshgrid(1:size(mask,1), 1:size(mask,2));
C       = [reshape(X,[],1), reshape(Y, [],1)];
mask_dist = [];
for i = 1:size(C,1)
    mask_dist(i) = min(sqrt(sum((C(reshape(logical(mask),[],1),:) - C(i,:))'.^2)));
end
%mask_dist = reshape(mask_dist, size(mask,1), size(mask,2));
mask_dist = mask_dist / max(mask_dist);

fig_supp = figure;
for p = 1:2
    subplot(2,2,p)
    hold on
    plot(reshape(Rmse_gllim(:,:,p) .* mask, 1,[]), reshape(CI(:,:,p) .* mask,1,[]), '.')
    %plot(reshape(Rmse_gllim(:,:,p) .* ~mask, 1,[]), reshape(CI(:,:,p) .* ~mask,1,[]), '.')
    plot([0 1.1*max(reshape(Rmse_gllim(:,:,p) .* mask, 1,[]))], [0 1.1*max(reshape(Rmse_gllim(:,:,p) .* mask, 1,[]))], 'k--')
    axis square
    xlabel('RMSE'); ylabel('Predicted RMSE')
    title('IN the dictionary')
    
    subplot(2,2,2+p)
    hold on
    tmp_ci = reshape(CI(:,:,p),[],1);
    tmp_err = reshape(Rmse_gllim(:,:,p), [],1);
    for i = 1:10:size(tmp_ci,1)
        plot(tmp_err(i), tmp_ci(i), '.', 'Color',[1 1 1]*mask_dist(i))
    end
    plot([0 1.1*max(reshape(Rmse_gllim(:,:,p) .* mask, 1,[]))], [0 1.1*max(reshape(Rmse_gllim(:,:,p) .* mask, 1,[]))], 'k--')
    axis square
    xlabel('RMSE'); ylabel('Predicted RMSE')
    title('OUT of the dico')
end


