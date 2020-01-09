function[]=inverseSIR(Ldatax,Ldatay,Ldatanoise, maskwave,logregul,image,image_mask)



%
% Ldatax are spectra from the learning database
%
% Ldatay are parameters from the learning database
%
% Ldatanoise are spectra from the learning database. A multi Gaussian noise
% has been added in order to estimate the regularization parameter.
%
% maskwave are the indexes of the considered wavelengths for the study 
%
% logregul is the logarithm of the values assumed by regularization parameter 
%
% image are the spectra we want to reverse
%
% image_mask is masking the spectra we do not want to reverse. For example,
% if we concentrate on the CO2 area, we will have to work on the area
% determined by image_mask.spectels(:,:,2)==1
%





% Importation of the data cubes (.lbl)


cd ~/dataGRSIR/
inter1=cub_info(Ldatax);
Ldatax=cub_read2(inter1);
clear inter1
inter2=cub_info(Ldatay);
Ldatay=cub_read2(inter2);
clear inter2
inter3=cub_info(Ldatanoise);
Ldatanoise=cub_read2(inter3);
clear inter3
inter4=cub_info(image);
image=cub_read2(inter4);
inter5=cub_info(image_mask);
image_mask=cub_read2(inter5);
clear inter5




% Since we gave values for the logarithm of the regularisation parameter,
% we calculate here the tested regularization values 

regul=10.^(logregul);
[pregul,nregul]=size(regul);


% finding indexes of the spectra from the CO2 area 

[ind1,ind2]=find(image_mask.spectels(:,:,2)==1);
ind=find(image_mask.spectels_flat(:,2)==1);

% Identity matrix for Tikhonov regularization

dimmask=size(maskwave);
Mid=eye(dimmask(1));

% working with parameter 2

for (k=1:nregul)
    
    % first we calculate the GRSIR axis
    % GRSIRp2 is a list depending on k, the kth tested regularization
    % parameter
    % GRSIRp2(k) contains:
    %                       - regul: the value of the regularization
    %                           parameter
    %                       - tikhonov: a matrix with the GRSIR axes.
    %                           tikhonov(:,1) is the first GRSIR axis
    %                       - projx the projected spectra on the first GRSIR
    %                           axis
    %                       - y the associated parameters values
    
GRSIRp2(k).res=regtikhonov(Ldatax.spectels_flat(:,maskwave),Ldatay.spectels_flat(:,2),regul(k),Mid);
    
    % Here, we calcultate the average of the reduced spectra in each slice 
    
EstimGRSIRp2(k).fn=meanslice(GRSIRp2(k).res.projx,GRSIRp2(k).res.y);

    % Linear piecewise interpolation and estimations of the noisy spectra
    % from the learning database
    
%EstimGRSIRp2(k).estimx=estiminterp1(EstimGRSIRp2(k).fn.moy,EstimGRSIRp2(k).fn.param,Ldatax.spectels_flat(:,maskwave)*GRSIRp2(k).res.tikhonov(:,1));


EstimGRSIRp2(k).estimxnoise=estiminterp1(EstimGRSIRp2(k).fn.moy,EstimGRSIRp2(k).fn.param,Ldatanoise.spectels_flat(:,maskwave)*GRSIRp2(k).res.tikhonov(:,1));

% NRMSE 
VcritGRSIRp2.test(k)=regul(k);


%VcritGRSIRp2.NRMSE(k)=calculNRMSE(Ldatay.spectels_flat(:,2),EstimGRSIRp2(k).estimx);
VcritGRSIRp2.NRMSEnoise(k)=calculNRMSE(Ldatay.spectels_flat(:,2),EstimGRSIRp2(k).estimxnoise);


end


% Choice of the regularization parameter by minimizing the NRMSE

indp2=find(VcritGRSIRp2.NRMSEnoise==min(VcritGRSIRp2.NRMSEnoise));

% Calculation of the SIRC criterion for the chosen regularization parameter

SIRCp2=SIRC(Ldatax.spectels_flat(:,maskwave),Ldatay.spectels_flat(:,2),GRSIRp2(indp2).res.tikhonov(:,1));

% Inversion of the image

InversGRSIRp2_flat=estiminterp1(EstimGRSIRp2(indp2).fn.moy,EstimGRSIRp2(indp2).fn.param,image.spectels_flat((image_mask.spectels_flat(:,2)==1),maskwave)*GRSIRp2(indp2).res.tikhonov(:,1));
InversGRSIRp2=sparse(ind1,ind2,InversGRSIRp2_flat);


% saving results

save('~/res/SIRCp2','SIRCp2')
save('~/res/GRSIRp2','GRSIRp2')
save('~/res/EstimGRSIRp2','EstimGRSIRp2')
save('~/res/InversGRSIRp2_flat','InversGRSIRp2_flat')
save('~/res/VcritGRSIRp2','VcritGRSIRp2')
save('~/res/InversGRSIRp2','InversGRSIRp2')
save('~/res/indp2','indp2')

clear SIRCp2
clear GRSIRp2
clear EstimGRSIRp2
clear VcritGRSIRp2
clear InversGRSIRp2
clear indp2

% same operations for the others parameters

for (k=1:nregul)
GRSIRp3(k).res=regtikhonov(Ldatax.spectels_flat(:,maskwave),Ldatay.spectels_flat(:,3),regul(k),Mid);
EstimGRSIRp3(k).fn=meanslice(GRSIRp3(k).res.projx,GRSIRp3(k).res.y);
EstimGRSIRp3(k).estimx=estiminterp1(EstimGRSIRp3(k).fn.moy,EstimGRSIRp3(k).fn.param,Ldatax.spectels_flat(:,maskwave)*GRSIRp3(k).res.tikhonov(:,1));
EstimGRSIRp3(k).estimxnoise=estiminterp1(EstimGRSIRp3(k).fn.moy,EstimGRSIRp3(k).fn.param,Ldatanoise.spectels_flat(:,maskwave)*GRSIRp3(k).res.tikhonov(:,1));
VcritGRSIRp3.test(k)=regul(k);
VcritGRSIRp3.NRMSE(k)=calculNRMSE(Ldatay.spectels_flat(:,3),EstimGRSIRp3(k).estimx);
VcritGRSIRp3.NRMSEnoise(k)=calculNRMSE(Ldatay.spectels_flat(:,3),EstimGRSIRp3(k).estimxnoise);
end
indp3=find(VcritGRSIRp3.NRMSEnoise==min(VcritGRSIRp3.NRMSEnoise));
SIRCp3=SIRC(Ldatax.spectels_flat(:,maskwave),Ldatay.spectels_flat(:,3),GRSIRp3(indp3).res.tikhonov(:,1))
InversGRSIRp3_flat=estiminterp1(EstimGRSIRp3(indp3).fn.moy,EstimGRSIRp3(indp3).fn.param,image.spectels_flat((image_mask.spectels_flat(:,2)==1),maskwave)*GRSIRp3(indp3).res.tikhonov(:,1));
InversGRSIRp3=sparse(ind1,ind2,InversGRSIRp3_flat);

save('~/res/SIRCp3','SIRCp3')
save('~/res/GRSIRp3','GRSIRp3')
save('~/res/EstimGRSIRp3','EstimGRSIRp3')
save('~/res/InversGRSIRp3_flat','InversGRSIRp3_flat')
save('~/res/VcritGRSIRp3','VcritGRSIRp3')
save('~/res/InversGRSIRp3','InversGRSIRp3')
save('~/res/indp3','indp3')



clear SIRCp3
clear GRSIRp3
clear EstimGRSIRp3
clear VcritGRSIRp3
clear InversGRSIRp3
clear indp3



for (k=1:nregul)
GRSIRp4(k).res=regtikhonov(Ldatax.spectels_flat(:,maskwave),Ldatay.spectels_flat(:,4),regul(k),Mid);
EstimGRSIRp4(k).fn=meanslice(GRSIRp4(k).res.projx,GRSIRp4(k).res.y);
EstimGRSIRp4(k).estimx=estiminterp1(EstimGRSIRp4(k).fn.moy,EstimGRSIRp4(k).fn.param,Ldatax.spectels_flat(:,maskwave)*GRSIRp4(k).res.tikhonov(:,1));
EstimGRSIRp4(k).estimxnoise=estiminterp1(EstimGRSIRp4(k).fn.moy,EstimGRSIRp4(k).fn.param,Ldatanoise.spectels_flat(:,maskwave)*GRSIRp4(k).res.tikhonov(:,1));
VcritGRSIRp4.test(k)=regul(k);
VcritGRSIRp4.NRMSE(k)=calculNRMSE(Ldatay.spectels_flat(:,4),EstimGRSIRp4(k).estimx);
VcritGRSIRp4.NRMSEnoise(k)=calculNRMSE(Ldatay.spectels_flat(:,4),EstimGRSIRp4(k).estimxnoise);
end
indp4=find(VcritGRSIRp4.NRMSEnoise==min(VcritGRSIRp4.NRMSEnoise));
SIRCp4=SIRC(Ldatax.spectels_flat(:,maskwave),Ldatay.spectels_flat(:,4),GRSIRp4(indp4).res.tikhonov(:,1))
InversGRSIRp4_flat=estiminterp1(EstimGRSIRp4(indp4).fn.moy,EstimGRSIRp4(indp4).fn.param,image.spectels_flat((image_mask.spectels_flat(:,2)==1),maskwave)*GRSIRp4(indp4).res.tikhonov(:,1));
InversGRSIRp4=sparse(ind1,ind2,InversGRSIRp4_flat);

save('~/res/SIRCp4','SIRCp4')
save('~/res/GRSIRp4','GRSIRp4')
save('~/res/EstimGRSIRp4','EstimGRSIRp4')
save('~/res/InversGRSIRp4_flat','InversGRSIRp4_flat')
save('~/res/VcritGRSIRp4','VcritGRSIRp4')
save('~/res/InversGRSIRp4','InversGRSIRp4')
save('~/res/indp4','indp4')



clear SIRCp4
clear GRSIRp4
clear EstimGRSIRp4
clear VcritGRSIRp4
clear InversGRSIRp4
clear indp4




for (k=1:nregul)
GRSIRp5(k).res=regtikhonov(Ldatax.spectels_flat(:,maskwave),Ldatay.spectels_flat(:,5),regul(k),Mid);
EstimGRSIRp5(k).fn=meanslice(GRSIRp5(k).res.projx,GRSIRp5(k).res.y);
EstimGRSIRp5(k).estimx=estiminterp1(EstimGRSIRp5(k).fn.moy,EstimGRSIRp5(k).fn.param,Ldatax.spectels_flat(:,maskwave)*GRSIRp5(k).res.tikhonov(:,1));
EstimGRSIRp5(k).estimxnoise=estiminterp1(EstimGRSIRp5(k).fn.moy,EstimGRSIRp5(k).fn.param,Ldatanoise.spectels_flat(:,maskwave)*GRSIRp5(k).res.tikhonov(:,1));
VcritGRSIRp5.test(k)=regul(k);
VcritGRSIRp5.NRMSE(k)=calculNRMSE(Ldatay.spectels_flat(:,5),EstimGRSIRp5(k).estimx);
VcritGRSIRp5.NRMSEnoise(k)=calculNRMSE(Ldatay.spectels_flat(:,5),EstimGRSIRp5(k).estimxnoise);
end
indp5=find(VcritGRSIRp5.NRMSEnoise==min(VcritGRSIRp5.NRMSEnoise));
SIRCp5=SIRC(Ldatax.spectels_flat(:,maskwave),Ldatay.spectels_flat(:,5),GRSIRp5(indp5).res.tikhonov(:,1))
InversGRSIRp5_flat=estiminterp1(EstimGRSIRp5(indp5).fn.moy,EstimGRSIRp5(indp5).fn.param,image.spectels_flat((image_mask.spectels_flat(:,2)==1),maskwave)*GRSIRp5(indp5).res.tikhonov(:,1));
InversGRSIRp5=sparse(ind1,ind2,InversGRSIRp5_flat);


save('~/res/SIRCp5','SIRCp5')
save('~/res/GRSIRp5','GRSIRp5')
save('~/res/EstimGRSIRp5','EstimGRSIRp5')
save('~/res/InversGRSIRp5_flat','InversGRSIRp5_flat')
save('~/res/VcritGRSIRp5','VcritGRSIRp5')
save('~/res/InversGRSIRp5','InversGRSIRp5')
save('~/res/indp5','indp5')


clear SIRCp5
clear GRSIRp5
clear EstimGRSIRp5
clear VcritGRSIRp5
clear InversGRSIRp5
clear indp5


for (k=1:nregul)
GRSIRp6(k).res=regtikhonov(Ldatax.spectels_flat(:,maskwave),Ldatay.spectels_flat(:,6),regul(k),Mid);
EstimGRSIRp6(k).fn=meanslice(GRSIRp6(k).res.projx,GRSIRp6(k).res.y);
EstimGRSIRp6(k).estimx=estiminterp1(EstimGRSIRp6(k).fn.moy,EstimGRSIRp6(k).fn.param,Ldatax.spectels_flat(:,maskwave)*GRSIRp6(k).res.tikhonov(:,1));
EstimGRSIRp6(k).estimxnoise=estiminterp1(EstimGRSIRp6(k).fn.moy,EstimGRSIRp6(k).fn.param,Ldatanoise.spectels_flat(:,maskwave)*GRSIRp6(k).res.tikhonov(:,1));
VcritGRSIRp6.test(k)=regul(k);
VcritGRSIRp6.NRMSE(k)=calculNRMSE(Ldatay.spectels_flat(:,6),EstimGRSIRp6(k).estimx);
VcritGRSIRp6.NRMSEnoise(k)=calculNRMSE(Ldatay.spectels_flat(:,6),EstimGRSIRp6(k).estimxnoise);
end
indp6=find(VcritGRSIRp6.NRMSEnoise==min(VcritGRSIRp6.NRMSEnoise));
SIRCp6=SIRC(Ldatax.spectels_flat(:,maskwave),Ldatay.spectels_flat(:,6),GRSIRp6(indp6).res.tikhonov(:,1))
InversGRSIRp6_flat=estiminterp1(EstimGRSIRp6(indp6).fn.moy,EstimGRSIRp6(indp6).fn.param,image.spectels_flat((image_mask.spectels_flat(:,2)==1),maskwave)*GRSIRp6(indp6).res.tikhonov(:,1));
InversGRSIRp6=sparse(ind1,ind2,InversGRSIRp6_flat);



save('~/res/SIRCp6','SIRCp6')
save('~/res/GRSIRp6','GRSIRp6')
save('~/res/EstimGRSIRp6','EstimGRSIRp6')
save('~/res/InversGRSIRp6_flat','InversGRSIRp6_flat')
save('~/res/VcritGRSIRp6','VcritGRSIRp6')
save('~/res/InversGRSIRp6','InversGRSIRp6')
save('~/res/indp6','indp6')


clear SIRCp6
clear GRSIRp6
clear EstimGRSIRp6
clear VcritGRSIRp6
clear InversGRSIRp6
clear indp6

load('~/res/InversGRSIRp2_flat')
load('~/res/InversGRSIRp3_flat')
load('~/res/InversGRSIRp4_flat')
load('~/res/InversGRSIRp5')
load('~/res/InversGRSIRp6')
load('~/res/InversGRSIRp4')
load('~/res/InversGRSIRp3')
load('~/res/InversGRSIRp2')



% The sum of the proportions hase to be one. The proportion of dust is
% deduced from the other proportions. If some proportions of dust are
% negative then the the proportion of dust in estimated by GRSIR and the
% proportion of CO2 is deduced from the proportion of the water and CO2



InversGRSIRp2_flat_b=1-InversGRSIRp3_flat-InversGRSIRp4_flat;
indneg=find(InversGRSIRp2_flat_b<0);
InversGRSIRp2_flat_b(indneg)=InversGRSIRp2_flat(indneg);
InversGRSIRp3_flat_b=InversGRSIRp3_flat;
InversGRSIRp3_flat_b(indneg)=1-InversGRSIRp2_flat(indneg)-InversGRSIRp4_flat(indneg);
InversGRSIRp4_flat_b=InversGRSIRp4_flat;

InversGRSIRp1=sparse(ind1,ind2,Ldatay.spectels_flat(1,1));
InversGRSIRp2_b=sparse(ind1,ind2,InversGRSIRp2_flat_b);
InversGRSIRp3_b=sparse(ind1,ind2,InversGRSIRp3_flat_b);
InversGRSIRp4_b=sparse(ind1,ind2,InversGRSIRp4_flat_b);

save('~/res/InversGRSIRp2_flat_b','InversGRSIRp2_flat_b')
save('~/res/InversGRSIRp3_flat_b','InversGRSIRp3_flat_b')
save('~/res/InversGRSIRp4_flat_b','InversGRSIRp4_flat_b')


% Reconstructing a cube (.lbl) with GRSIR estimations

[t1,t2,t3]=size(InversGRSIRp1);
recons=zeros(t1,t2,t3);
inter4.cubename='recons';
recons(:,:,1)=InversGRSIRp1;
recons(:,:,2)=InversGRSIRp2_b;
recons(:,:,3)=InversGRSIRp3_b;
recons(:,:,4)=InversGRSIRp4_b;
recons(:,:,5)=InversGRSIRp5;
recons(:,:,6)=InversGRSIRp6;
cd ~/res/
save('-mat','~/res/recons.mat','recons')
tata=cub_create(recons);
cub_write(inter4,tata)




