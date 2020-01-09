function gllim_planeto_test_image_mrf(images,image_masks,suffix,beta)
% Utiliser suffix = 'Real345_K10_nonorm'
% EMAIL expliquant les bons parameteres : 2/07/2013
path(path,'../mvrvm/');
PATH='/local_scratch/deleforg/Data/Planeto/Training/';

%load([PATH,'normalizer',suffix]);
%load([PATH,'theta_PLVM-0',suffix]);
load([PATH,'theta_PLVM-2',suffix]);
%load([PATH,'theta_JGMM',suffix]);
% load([PATH,'theta_SIR1',suffix]);
% load([PATH,'theta_SIR2',suffix]);
% load([PATH,'theta_RVM',suffix]);

maxiter = 20; % Max number of MRF iterations 

%Lt=size(theta_PLVM0.c,1);
Lt=3;

% Add white as the lowest value in the colormap
clims=[0.9942, 0.9997; 0.0003, 0.0030; 50, 450];
% This value should be set to pixels out of the mask so that they appear white
dlims=clims(:,1)-(clims(:,2)-clims(:,1))./(size(colormap,1)-1);
imlims=[dlims,clims(:,2)];
Xtests=cell(1,3);
for i=1:2 % i=1:3
    image_masks{i}=image_masks{i}((end-100):end,:,:);    
    images{i}=images{i}((end-100):end,:,:);      
    mask=image_masks{i}(:,:,2);
    image_masks{i}=image_masks{i}(sum(mask,2)>0,sum(mask,1)>0,:);    
    images{i}=images{i}(sum(mask,2)>0,sum(mask,1)>0,:);
    [I,J,D]=size(images{i});
    Xtests{i}=reshape(images{i},[I*J,D])';                                         
end                                             

% Build image neighborhood matrix:
fprintf(1,'Building Sparse Neighborhood Matrices\n');
Nei=cell(1,3);
for i=1:2 % i=1:3
    Nei{i} = sparse_neighbour_matrix(image_masks{i}(:,:,2));
end  

% %%%%% PLVM Lw=0 %%%%%%%%          
% colormap([[1,1,1];colormap]);
% method_name='PLVM-0';
% fprintf(1,'PLVM-0\n');
% for i=1:2 % i=1:3     
%     [I,J,~]=size(images{i});
%     image_mask=image_masks{i};
%     Xtest=Xtests{i};
%     Ytest_est=gllim_inverse_map(Xtest,theta_PLVM0,0);  
%     Ytest_est=Ytest_est(1:Lt,:);
% %     Ytest_est=bsxfun(@plus,bsxfun(@times,Ytest_est,normalizer.stdY),...
% %                            normalizer.meanY);    
%     for l=1:2 % l=1:Lt
%         subplot(2,2,(i-1)*2+l); % subplot(3,3,(i-1)*3+l); 
%         fprintf(1,[method_name,': min:%f 2.2 max:%f 2.2\n'],...
%                                   min(Ytest_est(l,:)),max(Ytest_est(l,:)));        
%         Ytest_est(l,Ytest_est(l,:)<clims(l,1))=clims(l,1);
%         Ytest_est(l,Ytest_est(l,:)>clims(l,2))=clims(l,2);        
% %         mask=sign(Ytest_est(l,1)).*log(image_mask(:,:,2))+1;
%         img=reshape(Ytest_est(l,:),[I,J]);
%         img(image_mask(:,:,2)==0)=dlims(l);
%         imagesc(img,imlims(l,:));
% %         imagesc(reshape(Ytest_est(l,:),[I,J]).*mask,clims(l,:));           
% %         fig=figure;clf(fig);
% %         mask_lin=reshape(image_mask(:,:,2),[I*J,1]);
% %         hist(Ytest_est(l,mask_lin==1),50);    
%     end
% end
% 
% %%%%% PLVM Lw=2 Beta=0 %%%%%%%%         
% fig=figure;clf(fig);
% colormap([[1,1,1];colormap]);
% method_name='PLVM-2';
% fprintf(1,'PLVM-2 beta=0\n');
% for i=1:2 % i=1:3     
%     [I,J,~]=size(images{i});
%     image_mask=image_masks{i};
%     Xtest=Xtests{i}; 
%     Ytest_est=gllim_inverse_map(Xtest,theta_PLVM2,0);  
%     Ytest_est=Ytest_est(1:Lt,:);
% %     Ytest_est=bsxfun(@plus,bsxfun(@times,Ytest_est,normalizer.stdY),...
% %                            normalizer.meanY);    
%     for l=1:2 % l=1:Lt
%         subplot(2,2,(i-1)*2+l); % subplot(3,3,(i-1)*3+l);  
%         fprintf(1,[method_name,': min:%f 2.2 max:%f 2.2\n'],...
%             min(Ytest_est(l,:)),max(Ytest_est(l,:)));           
%         Ytest_est(l,Ytest_est(l,:)<clims(l,1))=clims(l,1);
%         Ytest_est(l,Ytest_est(l,:)>clims(l,2))=clims(l,2);  
% %         mask=sign(Ytest_est(l,1)).*log(image_mask(:,:,2))+1;
%         img=reshape(Ytest_est(l,:),[I,J]);
%         img(image_mask(:,:,2)==0)=dlims(l);
%         imagesc(img,imlims(l,:));
% %         imagesc(reshape(Ytest_est(l,:),[I,J]).*mask,clims(l,:));           
% %         fig=figure;clf(fig);
% %         mask_lin=reshape(image_mask(:,:,2),[I*J,1]);
% %         hist(Ytest_est(l,mask_lin==1),50);    
%     end
% end

%%%%% PLVM Lw=2 %%%%%%%%         
fig=figure;clf(fig);
colormap([[1,1,1];colormap]);
method_name='PLVM-2';
fprintf(1,'PLVM-2 beta=%g\n',beta);
for i=1:2 % i=1:3     
    [I,J,~]=size(images{i});
    image_mask=image_masks{i};
    Xtest=Xtests{i};    
    [Ytest_est,alpha]=gllim_inverse_map_mrf(Xtest,theta_PLVM2,beta,Nei{i},reshape(image_mask(:,:,2),[I*J,1]),maxiter,2);
    [~,maxalpha] = max(alpha,[],2); % N*1
    img_alpha = reshape(maxalpha,[I,J]);
    img_alpha(~image_mask(:,:,2))=-1;
    subplot(2,3,(i-1)*3+1);
    imagesc(img_alpha,[-1,10]);
    Ytest_est=Ytest_est(1:Lt,:);
%     Ytest_est=bsxfun(@plus,bsxfun(@times,Ytest_est,normalizer.stdY),...
%                            normalizer.meanY);    
    for l=1:2 % l=1:Lt
        subplot(2,3,(i-1)*3+l+1); % subplot(3,3,(i-1)*3+l);  
        fprintf(1,[method_name,': min:%f 2.2 max:%f 2.2\n'],...
            min(Ytest_est(l,:)),max(Ytest_est(l,:)));           
        Ytest_est(l,Ytest_est(l,:)<clims(l,1))=clims(l,1)+1e-04;
        Ytest_est(l,Ytest_est(l,:)>clims(l,2))=clims(l,2);  
%         mask=sign(Ytest_est(l,1)).*log(image_mask(:,:,2))+1;
        img=reshape(Ytest_est(l,:),[I,J]);
        img(image_mask(:,:,2)==0)=dlims(l);
        imagesc(img,imlims(l,:));
%         imagesc(reshape(Ytest_est(l,:),[I,J]).*mask,clims(l,:));           
%         fig=figure;clf(fig);
%         mask_lin=reshape(image_mask(:,:,2),[I*J,1]);
%         hist(Ytest_est(l,mask_lin==1),50);    
    end
end

% %%%%% Joint GMM %%%%%%%%
% fig=figure;clf(fig);
% colormap([[1,1,1];colormap]);
% method_name='JGMM';
% for i=1:2 % i=1:3    
%     [I,J,~]=size(images{i});
%     image_mask=image_masks{i};
%     Xtest=Xtests{i};    
%     Ytest_est=gllim_inverse_map(Xtest,theta_JGMM,0);  
%     Ytest_est=Ytest_est(1:Lt,:);
% %     Ytest_est=bsxfun(@plus,bsxfun(@times,Ytest_est,normalizer.stdY),...
% %                            normalizer.meanY);    
%     for l=1:2 % l=1:Lt
%         subplot(2,2,(i-1)*2+l); % subplot(3,3,(i-1)*3+l); 
%         fprintf(1,[method_name,': min:%f 2.2 max:%f 2.2\n'],...
%             min(Ytest_est(l,:)),max(Ytest_est(l,:)));            
%         Ytest_est(l,Ytest_est(l,:)<clims(l,1))=clims(l,1);
%         Ytest_est(l,Ytest_est(l,:)>clims(l,2))=clims(l,2);          
% %         mask=sign(Ytest_est(l,1)).*log(image_mask(:,:,2))+1;
%         img=reshape(Ytest_est(l,:),[I,J]);
%         img(image_mask(:,:,2)==0)=dlims(l);
%         imagesc(img,imlims(l,:));
% %         imagesc(reshape(Ytest_est(l,:),[I,J]).*mask,clims(l,:));          
% %         fig=figure;clf(fig);
% %         mask_lin=reshape(image_mask(:,:,2),[I*J,1]);
% %         hist(Ytest_est(l,mask_lin==1),50);    
%     end
% end
% 
% for i=1:2 % i=1:3
%     [I,J,D]=size(images{i});
%     images{i}=bsxfun(@rdivide,bsxfun(@minus,images{i},reshape(normalizer.meanX,[1,1,D])),...
%                                                       reshape(normalizer.stdX,[1,1,D]));
%     Xtests{i}=reshape(images{i},[I*J,D])';                                                  
% end      

% %%% SIR Q=1 %%%
% fig=figure;clf(fig);
% colormap([[1,1,1];colormap]);
% method_name='SIR1';
% for i=1:2 % i=1:3  
%     [I,J,~]=size(images{i});
%     image_mask=image_masks{i};
%     Xtest=Xtests{i}; 
%     for l=1:2 % l=1:Lt
%         subplot(2,2,(i-1)*2+l); % subplot(3,3,(i-1)*3+l);          
%         betal=reshape(theta_SIR1.beta(:,1,l),[D,1]);
%         Ytest_est = feval(theta_SIR1.cfit,Xtest'*betal)';
%         Ytest_est=Ytest_est*normalizer.stdY(l)+normalizer.meanY(l);   
%         fprintf(1,[method_name,': min:%f 2.2 max:%f 2.2\n'],...
%             min(Ytest_est),max(Ytest_est));         
%         Ytest_est(Ytest_est<clims(l,1))=clims(l,1);
%         Ytest_est(Ytest_est>clims(l,2))=clims(l,2); 
% %         mask=sign(Ytest_est(1)).*log(image_mask(:,:,2))+1;
%         img=reshape(Ytest_est,[I,J]);
%         img(image_mask(:,:,2)==0)=dlims(l);
%         imagesc(img,imlims(l,:));
% %         imagesc(reshape(Ytest_est,[I,J]).*mask,clims(l,:));         
% %         fig=figure;clf(fig);
% %         mask_lin=reshape(image_mask(:,:,2),[I*J,1]);
% %         hist(Ytest_est(mask_lin==1),50);    
%     end
% end
% 
% %%% SIR Q=2 %%%
% fig=figure;clf(fig);
% colormap([[1,1,1];colormap]);
% method_name='SIR2';
% for i=1:2 % i=1:3   
%     [I,J,~]=size(images{i});
%     image_mask=image_masks{i};
%     Xtest=Xtests{i};     
%     for l=1:2 % l=1:Lt
%         subplot(2,2,(i-1)*2+l); % subplot(3,3,(i-1)*3+l);        
%         betal=reshape(theta_SIR2.beta(:,:,l),[D,2]);
%         Ytest_est = feval(theta_SIR2.cfit,Xtest'*betal(:,1),Xtest'*betal(:,2));
%         Ytest_est=Ytest_est*normalizer.stdY(l)+normalizer.meanY(l); 
%         fprintf(1,[method_name,': min:%f 2.2 max:%f 2.2\n'],...
%             min(Ytest_est),max(Ytest_est));        
%         Ytest_est(Ytest_est<clims(l,1))=clims(l,1);
%         Ytest_est(Ytest_est>clims(l,2))=clims(l,2);    
% %         mask=sign(Ytest_est(1)).*log(image_mask(:,:,2))+1;
%         img=reshape(Ytest_est,[I,J]);
%         img(image_mask(:,:,2)==0)=dlims(l);
%         imagesc(img,imlims(l,:));
% %         imagesc(reshape(Ytest_est,[I,J]).*mask,clims(l,:));             
% %         fig=figure;clf(fig);
% %         mask_lin=reshape(image_mask(:,:,2),[I*J,1]);
% %         hist(Ytest_est(mask_lin==1),50);   
%     end
% end
%        
% %%% RVM %%%
% fig=figure;clf(fig);
% colormap([[1,1,1];colormap]);
% method_name='RVM';
% for i=1:2 % i=1:3  
%     [I,J,~]=size(images{i});
%     image_mask=image_masks{i};
%     Xtest=Xtests{i};     
%     PHI	= sbl_kernelFunction(Xtest',theta_RVM.Xtrain',theta_RVM.kernel,theta_RVM.width);
%     Ytest_est	= (PHI*theta_RVM.weights)';  
%     Ytest_est = bsxfun(@plus,bsxfun(@times,Ytest_est,normalizer.stdY),...
%                            normalizer.meanY);
%     for l=1:2 % l=1:Lt
%         subplot(2,2,(i-1)*2+l); % subplot(3,3,(i-1)*3+l); 
%         fprintf(1,[method_name,': min:%f 2.2 max:%f 2.2\n'],...
%             min(Ytest_est(l,:)),max(Ytest_est(l,:)));          
%         Ytest_est(l,Ytest_est(l,:)<clims(l,1))=clims(l,1);
%         Ytest_est(l,Ytest_est(l,:)>clims(l,2))=clims(l,2); 
% %         mask=sign(Ytest_est(l,1)).*log(image_mask(:,:,2))+1;
%         img=reshape(Ytest_est(l,:),[I,J]);
%         img(image_mask(:,:,2)==0)=dlims(l);
%         imagesc(img,imlims(l,:));
% %         imagesc(reshape(Ytest_est(l,:),[I,J]).*mask,clims(l,:));             
% %         fig=figure;clf(fig);
% %         mask_lin=reshape(image_mask(:,:,2),[I*J,1]);
% %         hist(Ytest_est(l,mask_lin==1),50);  
%     end
% end

end