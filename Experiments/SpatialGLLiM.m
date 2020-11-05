
addpath(genpath('functions'))
addpath(genpath('tools'))

clear

%%

pb   	= 2;

params	= [1 2 3 4];
snr  	= 10;

N    	= 128;

nb_train_signals = 3000;

Parameters = [];
Parameters.K = 150;
Parameters.cstr.Sigma  = 'd*';
Parameters.cstr.Gammat = ''; 
Parameters.cstr.Gammaw = '';
Parameters.Lw = 0;
snr_train  	= 60;

int     = [.01 1];
p       = [.01 .01];
while min(abs(pdist(p',@(x,y) x-y))) < .05, p = 0.1 + 0.9*rand(1,10); end


%% Create test data

ph = phantom(N);
mask = zeros(size(ph));
u = unique(ph);
for i = 1:numel(u)
    mask(ph == u(i)) = i;
end

%%

val = rand(numel(u), numel(params));

for t = 1:numel(params)
    for i = 1:size(ph,1)
        for j = 1:size(ph,2)
            Yt(i,j,t) = val(mask(i,j),t);
        end
    end
end

Ytest = reshape(Yt, size(Yt,1) * size(Yt,2), size(Yt,3));
Xtest = [];
parfor i = 1:size(Ytest,1)
    Xtest(i,:) = toyMRsignal(Ytest(i,:), p(1:numel(params)));
end

artefact = ones(size(Yt,1), size(Yt,2));
switch pb
    case 1 %Simple bands
        artefact(20:22,:) = 0.25;
        artefact(31:35,:) = 0.25;
        artefact(40,:) = 0.25;
    case 2 % Sine bands
        artefact = 0.6 + 0.4*repmat(sin((1:size(ph,1))*8*pi/size(ph,1)), size(ph,1),1);
    case 3 %Salt and pepper noise
        artefact = imnoise(artefact,'salt & pepper');
        artefact(artefact == 0) = 0.1;
    case 4 % 
        Xtest = reshape(Xtest,N,N,size(Xtest,2));
        t = zeros(size(Xtest));
        t(end/2+1:end,:) = 0.5 * Xtest(1:end/2,:);
        Xtest = Xtest + t;
        Xtest = reshape(Xtest, N*N,size(Xtest,3));
end
SNR = snr * ones(size(Yt,1), size(Yt,2)) .* artefact';


%% Create train data

Ytrain  	= int(1) + (int(2) - int(1)) * net(scramble(sobolset(numel(params)),'MatousekAffineOwen'),nb_train_signals);
parfor sim = 1:size(Ytrain,1)
    Xtrain(sim,:) = toyMRsignal(Ytrain(sim,:),p(1:numel(params)));
end
DicoR{1}.MRSignals = AddNoise(abs(Xtrain), snr_train); 
DicoR{1}.Parameters.Par = Ytrain;


%% Train

[~,Parameters] = AnalyzeMRImages([], DicoR, 'DBL', Parameters);


%% Estimate

% If N the number of observations:
% - Nei (N x N)           % Sparse binary matrix, Nei(n,m)=1 iff y(:,n) and
%                           y(:,m) are neighboring observations
% - mask (N x 1)          % Binary, mask(n)=1 iff n is within the mask

verb    = 0;

m       = reshape(ones(size(ph)),1,[]); % TODO: find the correct slice
% patch   = [0 1 0;
%            1 1 1;
%            0 1 0];
% patch   = [1 zeros(1,size(ph,1)-2) 1 1 1 zeros(1,size(ph,1)-2) 1];
patch   = [1 1 1;
           1 1 1;
           1 1 1];
patch   = [1 1 1 zeros(1,size(ph,1)-3) 1 1 1 zeros(1,size(ph,1)-3) 1 1 1];

% patch   = [0 0 1 0 0;
%            0 1 1 1 0;
%            1 1 1 1 1;
%            0 1 1 1 0;
%            0 0 1 0 0;]
% patch   = [0 0 1 0 0 zeros(1,size(ph,1)-5) 0 1 1 1 0 zeros(1,size(ph,1)-5) 1 1 1 1 1 zeros(1,size(ph,1)-5) 0 1 1 1 0 zeros(1,size(ph,1)-5) 0 0 1 0 0;];

Nei = speye(size(Ytest,1));
parfor i = 1:size(Nei,1)
    Nei(i,:) = conv2(full(Nei(i,:)),patch,'same');
end

Xtest_noisy = AddNoise(Xtest, reshape(SNR,[],1));


theta       = Parameters.theta;

beta1       = 10;
Yspatial_b1 = gllim_inverse_map_mrf(Xtest_noisy', theta, beta1, Nei, m, Parameters.maxiter, verb);
Yspatial_b1	= (Yspatial_b1' .* Parameters.factors.Ystd) + Parameters.factors.Ymean;

beta2       = 50;
Yspatial_b2 = gllim_inverse_map_mrf(Xtest_noisy', theta, beta2, Nei, m, Parameters.maxiter, verb);
Yspatial_b2	= (Yspatial_b2'.* Parameters.factors.Ystd) + Parameters.factors.Ymean;

beta3       = 200;
Yspatial_b3 = gllim_inverse_map_mrf(Xtest_noisy', theta, beta3, Nei, m, Parameters.maxiter, verb);
Yspatial_b3	= (Yspatial_b3'.* Parameters.factors.Ystd) + Parameters.factors.Ymean;

Ygllim      = gllim_inverse_map_and_cov(Xtest_noisy', theta, verb);
Ygllim      = (Ygllim'.* Parameters.factors.Ystd) + Parameters.factors.Ymean;


%

bounds = {[0 1], [0 1]};
bounds_err = {[0 .5],[0 .5]};%{[0 .15]./3, [0 30e-6]./3};

for i = 1:2
    
    figure
    h(1) = subplot(2,5,1);
    imagesc(reshape(Ytest(:,i),size(ph,1),size(ph,2)), bounds{i});
    axis off; title('Original')
    
    %Ygllim
    h(2) = subplot(2,5,2);
    imagesc(reshape(Ygllim(:,i),size(ph,1),size(ph,2)), bounds{i});
    axis off; title('DBL')
    %Yspatial
    h(3) = subplot(2,5,3);
    imagesc(reshape(Yspatial_b1(:,i),size(ph,1),size(ph,2)), bounds{i});
    axis off; title(['DBL-Spatial - beta = ' num2str(beta1)])
    h(4) = subplot(2,5,4);
    imagesc(reshape(Yspatial_b2(:,i),size(ph,1),size(ph,2)), bounds{i});
    axis off; title(['DBL-Spatial - beta = ' num2str(beta2)])
    h(5) = subplot(2,5,5);
    imagesc(reshape(Yspatial_b3(:,i),size(ph,1),size(ph,2)), bounds{i});
    axis off; title(['DBL-Spatial - beta = ' num2str(beta3)])

    %AE gllim
    g(1) = subplot(2,5,7);
    ae  = abs(reshape(Ytest(:,i),size(ph,1),size(ph,2)) - reshape(Ygllim(:,i),size(ph,1),size(ph,2)));
    imagesc(ae, bounds_err{i});
    axis off; title('AError DBL')
    %AE spatial
    g(2) = subplot(2,5,8);
    imagesc(abs(reshape(Ytest(:,i),size(ph,1),size(ph,2)) - reshape(Yspatial_b1(:,i),size(ph,1),size(ph,2))), bounds_err{i});
    axis off; title('AError DBL-Spatial')
    g(3) = subplot(2,5,9);
    imagesc(abs(reshape(Ytest(:,i),size(ph,1),size(ph,2)) - reshape(Yspatial_b2(:,i),size(ph,1),size(ph,2))), bounds_err{i});
    axis off; title('AError DBL-Spatial')
    g(4) = subplot(2,5,10);
    imagesc(abs(reshape(Ytest(:,i),size(ph,1),size(ph,2)) - reshape(Yspatial_b3(:,i),size(ph,1),size(ph,2))), bounds_err{i});
    axis off; title('AError DBL-Spatial')
    
%     linkaxes(g,'z')    
%     linkaxes(h,'z')

    colormap(g(1),jet);
    colormap(g(2),jet);
    colormap(g(3),jet);
    colormap(g(4),jet);
end



%%

figure
beta = [0 5 10 50 100 200 500];

for b = 1:numel(beta)
    Yspatial_b1 = gllim_inverse_map_mrf(Xtest_noisy', theta, beta(b), Nei, m, Parameters.maxiter, verb);
    Yspatial_b1	= (Yspatial_b1' .* Parameters.factors.Ystd) + Parameters.factors.Ymean;
    Rmse(b,:) = EvaluateEstimation(Ytest,Yspatial_b1);
end

plot(beta,mean(Rmse'))