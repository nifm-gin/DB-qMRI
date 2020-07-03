function [x_exp,alpha] = gllim_inverse_map_mrf(y,theta,beta,Nei,mask,maxiter,verb)
%%%%%%%%%%%%%%%%% Inverse Mapping from Gllim Parameters %%%%%%%%%%%%%%%%%%%
%%%% Author: Antoine Deleforge (July 2012) - antoine.deleforge@inria.fr %%%
% Description: Map N observations y using the inverse conditional
% expectation E[x|y;theta] of the gllim model with parameters theta.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Input %%%%
% - y (DxN)               % Input observations to map
% - theta  (struct)       % Gllim model parameters
%   - theta.c (LxK)       % Gaussian means of X's prior
%   - theta.Gamma (LxLxK) % Gaussian covariances of X's prior
%   - theta.pi (1xK)      % Gaussian weights of X's prior
%   - theta.A (DxLxK)     % Affine transformation matrices
%   - theta.b (DxK)       % Affine transformation vectors
%   - theta.Sigma (DxDxK) % Error covariances
% - beta (>0)             % Smoothness parameters
% - Nei (N x N)           % Sparse binary matrix, Nei(n,m)=1 iff y(:,n) and
%                           y(:,m) are neighboring observations
% - mask (N x 1)          % Binary, mask(n)=1 iff n is within the mask
% - maxiter (int)         % Maximum number of iterations
% - verb {0,1,2}          % Verbosity (def 1)
%%%% Output %%%%
% - x_exp (LxN)           % Posterior mean estimates E[xn|yn;theta]
% - alpha (NxK)           % Weights of the posterior GMMs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[D,N]=size(y);
[L,K]=size(theta.c);
% ======================Inverse density parameters=========================
if(verb>=1);fprintf(1,'Compute c_k*, Gamma_k* and log(alpha_kn)\n');end;
cs = zeros(D,K);
Gammas = zeros(D,D,K);
logpyZ=zeros(N,K); % Log-Conditionals log(pik*p(y|Z=k;theta))  
for k=1:K
    if(verb>=2), fprintf(1,'  k=%d ',k); end
    if(verb>=2), fprintf(1,'AbcG '); end
    Ak=reshape(theta.A(:,:,k),D,L); % DxL
    bk=reshape(theta.b(:,k),D,1); % Dx1
    Sigmak=reshape(theta.Sigma(:,:,k),D,D); %DxD
    ck=reshape(theta.c(:,k),L,1); % Lx1
    Gammak=reshape(theta.Gamma(:,:,k),L,L); % LxL 
    
    if(verb>=2), fprintf(1,'cks '); end
    cs(:,k)=Ak*ck+bk; % Dx1
    
    if(verb>=2), fprintf(1,'Gks '); end
    Gammas(:,:,k)=Sigmak+Ak*Gammak*Ak'; % DxD
    
    if(verb>=2), fprintf(1,'logalpha \n'); end 
    logpyZ(:,k)=log(theta.pi(k))+...
                  loggausspdf(y,cs(:,k),Gammas(:,:,k)); % Nx1    
end
    

if(verb>=1);fprintf(1,'Compute log(qZn_prior(k))\n');end;
logqZ = repmat(log(theta.pi),[N,1]); % Log-Z prior estimate
logqZ(mask==0,:) = NaN;

[~,maxidx] = max(logqZ(mask==1,:),[],2);
prevmeanq = mean(maxidx);
converged = false;
iter=0;

display = 0;
I = 8; J = 8;
if display == 1
    [~,kcolor]=max(logpyZ,[],2);
    kcolor=exp(logpyZ(:,2));
    figure(42);clf(42);
    imagesc(reshape(kcolor,[I,J])) 
end


while(iter < maxiter && ~converged)
    iter=iter+1;
    if(verb>=2)
        fprintf(1,'iteration %d\n',iter);        
    end  
    
    % Asynchronuous method (slow, good convergence)
    p=randperm(N);
    for n=1:N
        if display == 1
           [~,kcolor]=max(logqZ,[],2);
           kcolor=exp(logqZ(:,2));            
           figure(43);clf(43);
           imagesc(reshape(kcolor,[I,J]));                        
         end
        if (n/1000==round(n/1000) && verb>=2)
            [~,maxidx] = max(logqZ(mask==1,:),[],2);
            fprintf(1,'    signal %d/%d, mean(logqZ)=%g min(logqZ)=%g\n',n,N,...
                      mean(maxidx),...
                      min(min(logqZ(isfinite(logqZ)))));            
        end        
        if mask(p(n))
            neighb = find(Nei(p(n),:));
            if ~isempty(neighb)
                logqZ(p(n),:) = logpyZ(p(n),:) + beta*logsumexp(logqZ(neighb,:));
            else
                logqZ(p(n),:) = logpyZ(p(n),:);
            end
            logqZ(p(n),logqZ(p(n),:)<-1e05)=-Inf;
%             check=logqZ(p(n),:);
            logqZ(p(n),:) = logqZ(p(n),:) - logsumexp(logqZ(p(n),:),2); % Normalization  
        end
    end 
    
%     % Synchronuous method (fast, poor convergence)
%     if(verb>=2);
%         fprintf(1,'  mean(mean(logqZ))=  %g\n',mean(mean(logqZ)));
%     end
%     [i,j,~] = find(Nei);
%     for k=1:K
%         splogqZk = sparse(i,j,logqZ(j,k),N,N);
%         logqZ(:,k) = logpyZ(:,k) + beta*splogsumexp(splogqZk,2);
%     end
%     logqZ = bsxfun(@minus,logqZ,logsumexp(logqZ,2)); % Normalization    

    [~,maxidx] = max(logqZ(mask==1,:),[],2);
    meanq=mean(maxidx);
    if(verb>=2)
        fprintf(1,'  delta=  %g\n',abs(prevmeanq-meanq)./abs(prevmeanq));
    end
    
    converged=abs(prevmeanq-meanq)./abs(prevmeanq)<=1;%1/sum(mask==1);
    prevmeanq=meanq;
end

% Computer the Z prior estimate:
logpZ = logqZ - logpyZ; % N x K
logpZ = bsxfun(@minus,logpZ,logsumexp(logpZ,2)); % Normalization  

if(verb>=1);fprintf(1,'Compute K projections to X space and weights\n');end;
proj=NaN(L,N,K);     % K projections to X space
logalpha=zeros(N,K); % Conditional log-weights log(p(Z=k|y;theta))   
for k=1:K
    if(verb>=2), fprintf(1,'  k=%d ',k); end
    if(verb>=2), fprintf(1,'AbcG '); end
    Ak=reshape(theta.A(:,:,k),D,L); % DxL
    bk=reshape(theta.b(:,k),D,1); % Dx1
    Sigmak=reshape(theta.Sigma(:,:,k),D,D); %DxD
    ck=reshape(theta.c(:,k),L,1); % Lx1
    Gammak=reshape(theta.Gamma(:,:,k),L,L); % LxL 
    
    if(verb>=2), fprintf(1,'iSks '); end
    invSigmaks2=eye(L)+Gammak*Ak'/Sigmak*Ak; %Sigmaks=invSigmaks2\Gammak

    if(verb>=2), fprintf(1,'Aks '); end
    Aks=invSigmaks2\Gammak*Ak'/Sigmak;     

    if(verb>=2), fprintf(1,'bks '); end
    bks=invSigmaks2\(ck-Gammak*Ak'/Sigmak*bk); 
    
    if(verb>=2), fprintf(1,'logalpha '); end 
    logalpha(:,k)=logpZ(:,k)'+...
                  loggausspdf(y,cs(:,k),Gammas(:,:,k)); % Nx1  
    
    if(verb>=2), fprintf(1,'projections '); end 
    proj(:,:,k)=bsxfun(@plus,Aks*y,bks); % LxN 

    if(verb>=2), fprintf(1,'\n'); end 
end
den=logsumexp(logalpha,2); % Nx1
logalpha=bsxfun(@minus,logalpha,den); % NxK Normalization
alpha=exp(logalpha); % NxK

x_exp=reshape(sum(bsxfun(@times,reshape(alpha,[1,N,K]),proj),3),L,N); %LxN

end
