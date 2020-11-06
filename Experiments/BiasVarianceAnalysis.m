


Init
disp(['Running experiment ' mfilename '.m'])


%clear
addpath(genpath('/home_eq/bouxfa/Code/Local/Experiments/Trajectories/functions'))
addpath(genpath('/home_eq/bouxfa/Code/Local/Experiments/Trajectories/inputs/vdspiral'))


%% Settings

display = 0;
nn_exec = 'cpu';
%load_model = [];

rng(4);

backup = 1;

snr = 40;

FOV     = 0.28;
res     = 128;

nb_param = 3;   % 2 for T1 and T2 estimates, and 3 for T1, T2 and Df estimates
nb_signals = [16^3]; %4096; 226981;


%% Create phantom

[ref, support, mask] = define_phantom(FOV, res);

if display == 1
    figure
    subplot(231)
    imagesc(ref)
    colormap jet; axis image; axis off
    title('Labels')
end

T1 = nan(size(ref)); T1_ = (1500:300:3000) *1e-3; %(ms)
T2 = nan(size(ref)); T2_ = (60:60:1000) *1e-3; %(ms)
df = nan(size(ref));

%relaxation times
u = unique(ref);
u(u == 0) = [];
for r = 1:numel(u)
    T1(ref == u(r)) = T1_(r);
    T2(ref == u(r)) = T2_(r);
end
T1 = T1 .* (1 + 1/10 * rand(size(T1)));
T2 = T2 .* (1 + 1/10 * rand(size(T2)));


%offresonance
max_df = 100 *pi/180; % (rad)
df_ = -max_df:(2*max_df)/size(df,1):max_df;
for i = 1:size(df,1)
    for j = 1:size(df,2)
        df(i,j) = df_(i) + df_(j);
    end
end

if display == 1
    subplot(234)
    imagesc(T1)
    colormap gray; axis image; axis off
    title('T1')
    subplot(235)
    imagesc(T2)
    colormap gray; axis image; axis off
    title('T2')
    subplot(236)
    imagesc(df)
    colormap gray; axis image; axis off
    title('df')
end


%% Simulate signals

% From Dan Ma's paper:

%simulation settings
FA      = (pi/180)* [90 ...
           10 + 50 *sin((2*pi/500) *(1:250)) + 5*randn(1,250) ...
           zeros(1,50) ...
           (10 + 50 *sin((2*pi/500) *(1:250)) + 5*randn(1,250)) /2 ...
           zeros(1,50) ...
           10 + 50 *sin((2*pi/500) *(1:250)) + 5*randn(1,250) ...
           zeros(1,50) ...
           (10 + 50 *sin((2*pi/500) *(1:99)) + 5*randn(1,99)) /2]; % Flip angles (rad)

TR      = perlin(length(FA)+100); TR = TR(1,:); % Repetition Time (sec)
TR      = filter(gausswin(50),1,TR);
TR      = TR(51:1050);
TR  =    1e-3 * (10.5 + 3.5* (TR - min(TR)) ./ max(TR - min(TR)));

S       = 1000;
sub     = 1;
TR      = TR(1:sub:S*sub);
FA      = FA(1:sub:S*sub);


if ~isempty(load_model)
    idx = 2;
    
    load(['outputs/' load_model], 'FA','TR', 'Params', 'NeuralNet','mx','mn','nb_signals');
    Params = Params{idx};
    NeuralNet = NeuralNet{idx};
    nb_signals = nb_signals(idx);
    mx = mx{idx};
    mn = mn{idx};
end

%img settings
N       = size(ref,1);

%simulation
if nb_param == 2
    Acq     = MRF_dictionary(reshape(T1,[],1), reshape(T2,[],1), [], FA, TR); 
elseif nb_param == 3
    Acq     = MRF_dictionary(reshape(T1,[],1), reshape(T2,[],1), reshape(df,[],1), FA, TR); 
end
Xacq    = (Acq.normalization.*Acq.magnetization).';
Iref    = reshape(Xacq, N,N, length(TR));

%add noise
Xacq(any(isnan(Xacq)'),:) = 0;

%generate images
Iacq    = reshape(Xacq, N,N, length(TR));

if display == 1
    figure
    subplot(211)
    plot(FA *180/pi); ylabel('FA (degrees)')
    subplot(212)
    plot(TR*1e3); ylabel('TR (ms)')
    xlabel('Repetitions')
end


%% Dico generation & Regression training

Parameters.K = 100;
Parameters.cstr.Sigma  = 'd*';
Parameters.cstr.Gammat = ''; 
Parameters.cstr.Gammaw = '';
Parameters.Lw = 0;
Parameters.maxiter = 250;
snr_train = 60;

int_T1  = 1e-3 * [200 3000];	% between 20 and 3000 ms
int_T2  = 1e-3 * [20  300];    % between 20 and 3000 ms
int_df  = pi/180 * [-200 200];  % +/- 400 Hz

    
clear X Y
nb_step = round(nb_signals^(1/nb_param));
v1  = int_T1(1) : (int_T1(2) - int_T1(1)) / (nb_step-1) : int_T1(2);
v2  = int_T2(1) : (int_T2(2) - int_T2(1)) / (nb_step-1) : int_T2(2);
v3  = int_df(1) : (int_df(2) - int_df(1)) / (nb_step-1) : int_df(2);
if nb_param == 2
    Y(:,1)  = repmat(v1, 1, length(v2));
    Y(:,2)  = repelem(v2, 1, length(v2));
    D       = MRF_dictionary(Y(:,1), Y(:,2), [], FA, TR); 
elseif nb_param == 3
    Y(:,1)  = repmat(v1, 1, length(v2)*length(v3));
    Y(:,2)  = repmat(repelem(v2, 1, length(v1)), 1,length(v3));
    Y(:,3)  = repelem(v3, 1, length(v2)*length(v3));
    D       = MRF_dictionary(Y(:,1), Y(:,2), Y(:,3), FA, TR); 
end
X       = (D.normalization.*D.magnetization).';
DicoG{1}.MRSignals = X; 
DicoG{1}.Parameters.Par = Y(:,1:nb_param);
disp('Dictionary for matching generated')


if isempty(load_model)
    % Compute training dataset
    clear X Y
    Y       = net(scramble(sobolset(nb_param),'MatousekAffineOwen'),nb_signals);
    Y(:,1)  = int_T1(1) + (int_T1(2) - int_T1(1)) * Y(:,1);
    Y(:,2)  = int_T2(1) + (int_T2(2) - int_T2(1)) * Y(:,2); 
    if nb_param == 2
        D       = MRF_dictionary(Y(:,1), Y(:,2), [], FA, TR); 
    elseif nb_param == 3
        Y(:,3)  = int_df(1) + (int_df(2) - int_df(1)) * Y(:,3);
        D       = MRF_dictionary(Y(:,1), Y(:,2), Y(:,3), FA, TR); 
    end 
    X       = (D.normalization.*D.magnetization).';
    X       = AddNoise(X, snr_train); 
    DicoR{1}.MRSignals = [real(X) imag(X)];
    DicoR{1}.Parameters.Par = Y(:,1:nb_param);

    if size(DicoG{1}.MRSignals,1) ~= size(DicoR{1}.MRSignals,1)
        warning('Sizes are not equals')
    end
    disp('Dictionary for learning generated')


    % Learn models
    [~, Params] = AnalyzeMRImages([],DicoR,'DBL',Parameters);
    disp('DB-SL model learnt')


    mn      = min(DicoR{1}.Parameters.Par );
    mx      = max(DicoR{1}.Parameters.Par  - mn);
    YtrainNN = (DicoR{1}.Parameters.Par  - mn) ./ mx;
    NeuralNet = EstimateNNmodel(DicoR{1}.MRSignals,YtrainNN,0,nn_exec,100,6);
    disp('DB-DL model learnt')
end


%%

% Generate test data
if nb_param == 2
    Ytest 	= [reshape(T1,[],1), reshape(T2,[],1)];
elseif nb_param == 3
    Ytest 	= [reshape(T1,[],1), reshape(T2,[],1), reshape(df,[],1)];
end
Xtest   = reshape(Iacq,[],size(Iacq,3));

nan_val = any(isnan(Ytest)');
Ytest(nan_val,:) = nan;
Xtest(nan_val,:) = nan;

mc_sim = 100;
clear *grid *gllim *nn
for mc = 1:mc_sim
    
    disp(mc)
    
    XtestN  = AddNoise(Xtest, snr);
    
    % Perform DBM
    Estim   = AnalyzeMRImages(XtestN,DicoG,'DBM',[],Ytest(:,1:nb_param));
    Ygrid(:,:,mc) = squeeze(Estim.GridSearch.Y);
    Rmse_grid(:,mc) = Estim.GridSearch.Errors.Rmse;

    % Perform DBSL
    Estim   = AnalyzeMRImages([real(XtestN) imag(XtestN)],[],'DBL', Params, Ytest(:,1:nb_param),[],snr);
%     Estim   = AnalyzeMRImages(real(XtestN),[],'DBL', Params, Ytest(:,1:nb_param),[],snr);
    Ygllim(:,:,mc) = squeeze(Estim.Regression.Y);
    Rmse_gllim(:,mc) = Estim.Regression.Errors.Rmse;
    CI_gllim(:,:,mc) = squeeze(Estim.Regression.Cov);
    
    % Perform DBDL
    Y     = EstimateParametersFromNNmodel([real(XtestN) imag(XtestN)], NeuralNet, nn_exec);
%     Y     = EstimateParametersFromNNmodel(real(XtestN), NeuralNet, nn_exec);
    Ynn(:,:,mc) = (Y .* mx) + mn;
    Rmse_nn(:,mc) = EvaluateEstimation(Ytest, Ynn(:,:,mc));
end

%%

for i = 1:size(Ytest,1)
    [Bias_grid(i,:), Var_grid(i,:)] = BiasVariance(Ytest(i,:), squeeze(Ygrid(i,:,:))');
    [Bias_gllim(i,:), Var_gllim(i,:)] = BiasVariance(Ytest(i,:), squeeze(Ygllim(i,:,:))');
    [Bias_nn(i,:), Var_nn(i,:)] = BiasVariance(Ytest(i,:), squeeze(Ynn(i,:,:))');
end
 

%%

if backup == 1
    save(['temp/' 'temp_bias_variance_analysis'])
end


%%

bds = {[0 .15], [0 .015], [0 50*pi/180]};

fig = figure;
plot_dbm = 1;
    
tmp = nanmean(Bias_gllim,3);
for i =1:size(Bias_gllim,1)
        if abs(Bias_gllim(i,1)) > 0.3
            Bias_gllim(i,1) = nan;
        end
        if abs(Bias_gllim(i,2)) > 0.015
            Bias_gllim(i,2) = nan;
        end
        if abs(Bias_gllim(i,3)) > 0.8
            Bias_gllim(i,3) = nan;
        end
end

tmp = nanmean(CI_gllim,3);

for p = 1:nb_param
    subplot(9,3,1+9*(p-1) +3)
    imagesc(reshape(abs(Bias_nn(:,p)),N,N), bds{p})
    axis image; colormap hot
    ylabel('DB-DL')
    subplot(9,3,2+9*(p-1) +3)
    imagesc(reshape(Var_nn(:,p),N,N).^.5, bds{p})
    axis image; axis off; colormap hot
    subplot(9,3,3+9*(p-1) +3)
    imagesc((reshape(Bias_nn(:,p).^2 + Var_nn(:,p),N,N)).^.5, bds{p})
    axis image; axis off; colormap hot
    

    subplot(9,3,1+9*(p-1) +6)
    imagesc(reshape(abs(Bias_gllim(:,p)),N,N), bds{p})
    axis image; colormap hot
    ylabel('DB-SL')
    subplot(9,3,2+9*(p-1) +6)
    imagesc(reshape(Var_gllim(:,p),N,N).^.5, bds{p})
    axis image; axis off; colormap hot; %colorbar
    subplot(9,3,3+9*(p-1) +6)
    imagesc((reshape(Bias_gllim(:,p).^2 + Var_gllim(:,p),N,N)).^.5, bds{p})
    axis image; axis off; colormap hot
    
    
    
    
    axis image; axis off; colormap hot
    if plot_dbm == 1
        subplot(9,3,1+9*(p-1))
        imagesc(reshape(abs(Bias_grid(:,p)),N,N), bds{p})
        axis image; colormap hot
        if p == 1, title('Bias'); end
        ylabel('DBM')
        subplot(9,3,2+9*(p-1))
        imagesc(reshape(Var_grid(:,p),N,N).^.5, bds{p})
        axis image; axis off; colormap hot;
        if p == 1, title('Var1/2'); end
        subplot(9,3,3+9*(p-1))
        imagesc((reshape(Bias_grid(:,p).^2 + Var_grid(:,p),N,N)).^.5, bds{p})
        axis image; axis off; colormap hot
        if p == 1, title('RMSE'); end
        
    end
    
    %table of mean result
    T(1,1,p) = nanmean(reshape(abs(Bias_grid(:,p)),1,[]));
    T(2,1,p) = nanmean(reshape(Var_grid(:,p).^.5,1,[]));
    T(3,1,p) = nanmean(reshape((Bias_grid(:,p).^2 + Var_grid(:,p)).^.5,1,[]));
    
    T(1,2,p) = nanmean(reshape(abs(Bias_nn(:,p)),1,[]));
    T(2,2,p) = nanmean(reshape(Var_nn(:,p).^.5,1,[]));
    T(3,2,p) = nanmean(reshape((Bias_nn(:,p).^2 + Var_nn(:,p)).^.5,1,[]));
    
    T(1,3,p) = nanmean(reshape(abs(Bias_gllim(:,p)),1,[]));
    T(2,3,p) = nanmean(reshape(Var_gllim(:,p).^.5,1,[]));
    T(3,3,p) = nanmean(reshape((Bias_gllim(:,p).^2 + Var_gllim(:,p)).^.5,1,[]));
end


T(:,:,1) = T(:,:,1) * 1e3;
T(:,:,2) = T(:,:,2) * 1e3;
T(:,:,3) = T(:,:,3) * 180/pi;


%%

if backup == 1
    savefig(fig, ['figures/' 'temp_bias_variance_analysis'])
end