%% m = MRDataAnalytical(phantom, s, w)
%
% Function that returns the MR data corresponding to the given phantom that
% is modulated by a sensitivity profile, for the given kspace samples.
%
% Matthieu Guerquin-Kern, Biomedical Imaging Group / EPF Lausanne,
% 31-10-2009 (dd-mm-yyyy)

function m = MRDataAnalytical(phantom,s,w)

k = diag(phantom.FOV)*w.'/2/pi;
m = zeros(1,size(k,2));
ind_ell = [];
ind_bez = [];
Nbez = [];
for p = 1:length(phantom.region)
    region = phantom.region{p};
    switch region.type
        case {'ellipse'}
            ind_ell = [ind_ell, p];
        case {'polygon'}
            ind_bez = [ind_bez, p];
            Nbez = [Nbez,size(region.vertex,1)];
        case {'bezier'}
            ind_bez = [ind_bez, p];
            Nbez = [Nbez,size(region.control,1)];
        otherwise
            error('%s: this kind of region is unknown',region.type);
    end
end
Nell = numel(ind_ell);
[Nbez,indsort] = sort(Nbez,'ascend');
ind_bez = ind_bez(indsort);clear indsort;

hbar = waitbar(0,'MR-DATA.m -> Estimating time...');
tcum = 0;
tell = 0;
tbez = 0;
if Nell>0
    t = clock();
    m = m + MRDataEllipse(phantom.region{ind_ell(1)},s,k);
    tell = etime(clock(),t);
    tcum = tcum + tell;
end
if numel(Nbez)>0
    t = clock();
    m = m + MRDataBezier(phantom.region{ind_bez(1)},s,k);
    t = etime(clock(),t);
    tbez = t/Nbez(1);
    tcum = tcum + t;
end
t_est = max(Nell*tell + sum(Nbez)*tbez,tcum);
waitbar(tcum/t_est,hbar,sprintf('MR-DATA.m -> Estimated remaning time: %.0f s (total time: %.0f)',t_est-tcum,t_est));
for i = 2:(numel(ind_ell))
    p = ind_ell(i);
    region = phantom.region{p};
    t = clock();
    m = m + MRDataEllipse(region,s,k);
    t = etime(clock(),t);
    tell = (i*tell+t)/(i+1);
    tcum = tcum+t;
    t_est = max(Nell*mean(tell) + sum(Nbez)*mean(tbez),tcum);
    waitbar(tcum/t_est,hbar,sprintf('MR-DATA.m -> Estimated remaning time: %.0f s (total time: %.0f)',t_est-tcum,t_est));
end
for i = 1:(numel(ind_bez)-1)
    p = ind_bez(end-i+1);
    region = phantom.region{p};
    t = clock();        
    m = m + MRDataBezier(region,s,k);
    t = etime(clock(),t);
    tbez = (i*tbez+t/Nbez(end-i+1))/(i+1);
    tcum = tcum+t;
    t_est = max(Nell*mean(tell) + sum(Nbez)*tbez,tcum);
    waitbar(tcum/t_est,hbar,sprintf('MR-DATA.m -> Estimated remaning time: %.0f s (total time: %.0f)',t_est-tcum,t_est));
end
close(hbar);
m = (exp(1i*pi*(k(1,:)+k(2,:))).*m).';
