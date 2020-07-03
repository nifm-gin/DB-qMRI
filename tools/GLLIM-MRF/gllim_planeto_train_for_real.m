function gllim_planeto_train_for_real(Xtrain,Ytrain,K,cstrS,suffix)
PATH='/HOMES/deleforge/Data/Planeto/Training/';
path(path,'../');
path(path,'../mvrvm');

% Number and dimension
[Lt,~]=size(Ytrain);
[D,~]=size(Xtrain);

% PLVM constraints
cstr.Sigma=cstrS;  
%cstr.Gamma='i*';
maxiter=75;

% MV RVM params:
kernel = 'poly3';
width = 6;
maxIts	= 15;

% % Normalize data
% normalizer.meanY=mean(Ytrain,2);
% normalizer.stdY=std(Ytrain,1,2);
% normalizer.meanX=mean(Xtrain,2);
% normalizer.stdX=std(Xtrain,1,2);

% Normalize data
normalizer.meanY=zeros(Lt,1);
normalizer.stdY=ones(Lt,1);
normalizer.meanX=zeros(D,1);
normalizer.stdX=ones(D,1);


save([PATH,'normalizer',suffix],'normalizer');

Ytrain=bsxfun(@rdivide,bsxfun(@minus,Ytrain,normalizer.meanY),normalizer.stdY);
Xtrain=bsxfun(@rdivide,bsxfun(@minus,Xtrain,normalizer.meanX),normalizer.stdX);


% %%% RVM %%%
% fprintf(1,['==== RVM ',kernel,' (%g) ====  '],width);
% kernel_	= strcat('+',kernel); 
% PHI	= sbl_kernelFunction(Xtrain',Xtrain',kernel_,width); % Compute kernel
% [weights, used, ~, ~] = mvrvm(PHI,Ytrain',maxIts); % Run mvrvm
% used	= used - 1;
% if used(1)~=0; kernel_(1) = []; else used(1) = [];end;
% theta_RVM.weights=weights;
% theta_RVM.Xtrain=Xtrain(:,used);
% theta_RVM.kernel=kernel_;
% theta_RVM.width=width;
% save([PATH,'theta_RVM',suffix],'theta_RVM');        

% %%%% PLVM Lw=0 %%%%%%%%
% fprintf(1,'==== PLVM Lw=0 ====  \n');
% [theta_PLVM0,~]=gllim(Ytrain,Xtrain,K,'Lw',0,'cstr',cstr,'maxiter',maxiter,'verb',2);             
% save([PATH,'theta_PLVM-0',suffix],'theta_PLVM0');

%%%% PLVM Lw=2 %%%%%%%%
fprintf(1,'==== PLVM Lw=2 ====  \n');
[theta_PLVM2,~]=gllim(Ytrain,Xtrain,K,'Lw',2,'cstr',cstr,'maxiter',maxiter,'verb',2);   
save([PATH,'theta_PLVM-2',suffix],'theta_PLVM2');

% %%% Joint GMM %%%
% fprintf(1,'==== JGMM ====  \n'); 
% [~, model, ~, ~] = emgm([Ytrain;Xtrain],K,maxiter,1);
% theta_JGMM = gmmj2gllim(model,Lt);
% save([PATH,'theta_JGMM',suffix],'theta_JGMM');

% %%% SIR1 %%%
% fprintf(1,'==== SIR1 ====  \n');
% nslices=20;
% Q=1;
% theta_SIR1.beta=zeros(D,Q,Lt); % Q-dim principal subspaces of each parameter
% for l=1:Lt
%     [passage,lambda,~]=sir2(Xtrain',Ytrain(l,:)',nslices); % DxD
%     [~,p]=sort(diag(lambda),'descend');
%     theta_SIR1.beta(:,:,l)=passage(:,p(1:Q)); % DxQ principal subspace
%     betal=reshape(theta_SIR1.beta(:,1,l),[D,1]);
%     theta_SIR1.cfit = fit(Xtrain'*betal,Ytrain(l,:)','poly3');
% end
% save([PATH,'theta_SIR1',suffix],'theta_SIR1');
% 
% %%% SIR2 %%%
% fprintf(1,'==== SIR2 ====  \n');
% nslices=20;
% Q=2;
% theta_SIR2.beta=zeros(D,Q,Lt); % Q-dim principal subspaces of each parameter
% for l=1:Lt
%     [passage,lambda,~]=sir2(Xtrain',Ytrain(l,:)',nslices); % DxD
%     [~,p]=sort(diag(lambda),'descend');
%     theta_SIR2.beta(:,:,l)=passage(:,p(1:Q)); % DxQ principal subspace
%     betal=reshape(theta_SIR2.beta(:,:,l),[D,Q]);
%     theta_SIR2.cfit = fit(Xtrain'*betal,Ytrain(l,:)','poly33');
% end
% save([PATH,'theta_SIR2',suffix],'theta_SIR2');

end