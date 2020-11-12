function [Yestim] = EstimateParametersFromAnalytic(t,tse, Xpre,Xpost, Params)

%We work in ms, thus convert sec to ms
t   = t * 1e3;
tse = tse * 1e3;


echotimes_used = (t <= tse/2);

for i = 1:size(Xpre,1)
    tmp = levenbergmarquardt('AB_t2s', t(echotimes_used),Xpre(i,echotimes_used),  [30 max(Xpre(i,echotimes_used))]);
    t2s_pre(i) = tmp(1);
    tmp = levenbergmarquardt('AB_t2s', t(echotimes_used),Xpost(i,echotimes_used), [30 max(Xpost(i,echotimes_used))]);
    t2s_post(i) = tmp(1);
end

echotimes_used = find(abs(t - tse) == min(abs(t - tse)));

DeltaR2star = 1 ./ t2s_post - 1 ./ t2s_pre;
DeltaR2     = log( Xpre(:,echotimes_used) ./ Xpost(:,echotimes_used) )' ./ tse;

if isempty(Params)
    gamma   = 2.675 *1e8;
    B0      = 4.7;
    DeltaKhi = 0.2785 *1e-6; %CGS
    ADC     = 800 *1e-12;
end

%BVf
Yestim(:,1) = ( 3 ./ (4*pi*gamma*B0*DeltaKhi) ) .* DeltaR2star;
%VSI %TODO: chek if the imaginary part can be extracted like this
Yestim(:,2) = 0.425 .* ( ADC ./ (gamma*B0*DeltaKhi) ).^.5 .* abs((DeltaR2star ./ DeltaR2).^(3/2));


