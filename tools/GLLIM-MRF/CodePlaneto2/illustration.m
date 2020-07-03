load XBaseAleatoireBruit.mat
load YBaseAleatoire.mat

x=XBaseAleatoireBruit;
y=YBaseAleatoire(:,3);

nbslices=20;
expoparamess=-10:1:0;
omega=cov(x);

nbess=length(expoparamess);
erreur=zeros(nbess,1);

for j=1:nbess;
	j
	param=10^(expoparamess(j));
	erreur(j)=cross_validation(x,y,10,param,nbslices,omega);
	erreur(j)
end;
figure;
plot(erreur);

[tmp, jopt]=min(erreur);
param=10^(expoparamess(jopt));

param=5*10^-7;
beta=myGRSIR(x,y,param,nbslices,omega);
xproj=x*beta;

figure
plot(xproj,y,'.');
