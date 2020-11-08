%% MRDataEllipse.m
%
% Returns the MR data corresponding to the indicator function of an ellipse
% that is modulated by a polynomial or sinusoidal profile, for the given
% kspace samples.
%
% Matthieu Guerquin-Kern, Biomedical Imaging Group / EPF Lausanne,
% 31-10-2009 (dd-mm-yyyy)

function m = MRDataEllipse(ellipse,s,k)

theta = ellipse.angle;
width = ellipse.width;
center = ellipse.center;
weight = ellipse.weight;

Rmat = [cos(theta) ,-sin(theta);sin(theta) ,cos(theta)];
Dmat = diag(width/2);

switch lower(s.model)
    case 'polynomial'
        m = weight*MRDataEllipsePolynomial(Dmat,Rmat,center,s.data,k);
    case 'sinusoidal'
        m = weight*MRDataEllipseSinusoidal(Dmat,Rmat,center,s.data,k);
    otherwise
        error('unknown sensitivity model. Must be either ''polynomial'' or ''sinusoidal''.')
end
end

%%
function m = MRDataEllipseSinusoidal(Dmat,Rmat,center,coeff,k)
N = length(coeff);
Nk = size(k,2);
L = floor(sqrt(N));

[kx,ky] = GenerateFullCart2DGrid(L*[1,1]);
kx = kx(:)/2;
ky = ky(:)/2;
kx = repmat(kx,1,Nk)+repmat(k(1,:),N,1);
ky = repmat(ky,1,Nk)+repmat(k(2,:),N,1);
wu = -2*pi*Dmat*Rmat'*[kx(:).';ky(:).'];
wux = reshape(wu(1,:),N,Nk);
wuy = reshape(wu(2,:),N,Nk);
modw = sqrt(wux.^2+wuy.^2);
Gval = G(1,modw);
m = 2*pi*det(Dmat)*coeff(:).'*(exp(2*pi*1j*(center(1)*kx+center(2)*ky)).*Gval);
end

%%
function m = MRDataEllipsePolynomial(Dmat,Rmat,center,coeff,k)
N = length(coeff);
D = floor((sqrt(1+N*8)-3)/2); % works only in 2D!

Np = N; % N is sufficient for matrix inversion. However, for high degree one should consider to increase this number in order to improve conditionning
u = zeros(2,Np);
u(1,:) = rand(1,Np)-0.5;
u(2,:) = rand(1,Np)-0.5;
r = Rmat*Dmat*u+repmat(center(:),[1,Np]);%figure(1);plot(r(1,:),r(2,:),'g.');hold on;plot(u(1,:),u(2,:),'r.');hold off;
Mr = Polynomial2DMatrix(r,D);
Mu = Polynomial2DMatrix(u,D);%cond(mu)
t = Mu\Mr*coeff(:);

lookup = BuildLookupG(D); % Build the lookup table for partial derivatives
m = zeros(1,size(k,2));
w = -2*pi*k;
wu = Dmat*Rmat'*w;
modw = sqrt(sum(wu.^2,1));
G_table = cell(1,D+1);
for i = 1:D+1
    G_table{i} = G(i,modw);
end

or = zeros(2,N);
for d = 0:D
    n = d*(d+1)/2;
    for i = 0:d
        or(:,n+i+1) = [i;d-i];
    end
end

for i = 1:N
    order = sum(or(:,i));
    table = lookup{order+1}{or(2,i)+1};
    for j = 1:length(table.n)
        m = m + (1j)^order*t(i)*G_table{table.n(j)}.*table.coeff(j).*wu(1,:).^table.o(1,j).*wu(2,:).^table.o(2,j);
    end
end
m = 2*pi*det(Dmat)*exp(-1j*center(:)'*w).*m;
end

%% BuildLookupG.m
%
% Build the lookup table for partial derivatives of G.
%
% SEE: G.m
%
% Matthieu Guerquin-Kern, Biomedical Imaging Group / EPF Lausanne,
% 31-10-2009 (dd-mm-yyyy)

function lookup = BuildLookupG(D)

if exist('lookupG.mat','file')
    load lookupG.mat;
else
    lookup = [];
end

if length(lookup)>=D+1
    tmp = lookup;
    lookup = cell(1,D+1);
    for d=1:D+1
        lookup{d} = tmp{d};
    end
    return;
end

if isempty(lookup)
    lookup{1}{1}.n = 1;
    lookup{1}{1}.o = [0;0];
    lookup{1}{1}.coeff = 1;
end

D_cur = length(lookup);
tmp = lookup;
lookup = cell(1,D+1);
for d=1:D_cur
    lookup{d} = tmp{d};
end

for d = D_cur:D
    lookup{d+1} = cell(1,d+1); % order of derivative wrt 2nd variable
    % the order wrt the second variable follows
    for d2 = 0:floor(d/2) % order of derivative wrt 2nd variable
        lookup{d+1}{d2+1}.n = [];
        lookup{d+1}{d2+1}.o = [];
        lookup{d+1}{d2+1}.coeff = [];
        for n = 1:length(lookup{d}{d2+1}.n)
            if lookup{d}{d2+1}.o(1,n) > 0
                lookup{d+1}{d2+1}.n     = [lookup{d+1}{d2+1}.n,     lookup{d}{d2+1}.n(n)];
                lookup{d+1}{d2+1}.o     = [lookup{d+1}{d2+1}.o,     lookup{d}{d2+1}.o(:,n)-[1;0]];
                lookup{d+1}{d2+1}.coeff = [lookup{d+1}{d2+1}.coeff, lookup{d}{d2+1}.coeff(n)*lookup{d}{d2+1}.o(1,n)];
            end
            lookup{d+1}{d2+1}.n     = [lookup{d+1}{d2+1}.n,     lookup{d}{d2+1}.n(n) + 1        ];
            lookup{d+1}{d2+1}.o     = [lookup{d+1}{d2+1}.o,     lookup{d}{d2+1}.o(:,n) + [1;0]  ];
            lookup{d+1}{d2+1}.coeff = [lookup{d+1}{d2+1}.coeff, -lookup{d}{d2+1}.coeff(n)       ];
        end
        % looking for simplifications
        for n = 1:length(lookup{d+1}{d2+1}.n)
            nn = lookup{d+1}{d2+1}.n(n);
            o1n = lookup{d+1}{d2+1}.o(1,n);
            ind = find(lookup{d+1}{d2+1}.n==nn&lookup{d+1}{d2+1}.o(1,:)==o1n);
            if numel(ind)>1
                if ind(1)==n
                    lookup{d+1}{d2+1}.coeff(n) = sum(lookup{d+1}{d2+1}.coeff(ind));
                    lookup{d+1}{d2+1}.coeff(ind(2:end)) = 0;
                end
            end
        end
        ind = find(lookup{d+1}{d2+1}.coeff~=0);
        [lookup{d+1}{d2+1}.n,ind2] = sort(lookup{d+1}{d2+1}.n(ind));
        lookup{d+1}{d2+1}.o = lookup{d+1}{d2+1}.o(:,ind(ind2));
        lookup{d+1}{d2+1}.coeff = lookup{d+1}{d2+1}.coeff(ind(ind2));
        
        % symmetry
        if d2~=d/2
            lookup{d+1}{d+1-d2}.n       =   lookup{d+1}{d2+1}.n;
            lookup{d+1}{d+1-d2}.o       =   zeros(size(lookup{d+1}{d2+1}.o));
            lookup{d+1}{d+1-d2}.o(1,:)  =   lookup{d+1}{d2+1}.o(2,:);
            lookup{d+1}{d+1-d2}.o(2,:)  =   lookup{d+1}{d2+1}.o(1,:);
            lookup{d+1}{d+1-d2}.coeff   =   lookup{d+1}{d2+1}.coeff;
        end
    end
end

save([pwd '/lookupG.mat'],'lookup');
end

%% G.m
%
% Function that computes the MRI response of a centered circle.
%
% Matthieu Guerquin-Kern, Biomedical Imaging Group / EPF Lausanne,
% 19-09-2009 (dd-mm-yyyy)

function val = G(n,z)
warning('OFF','all');
ind_small   = find(abs(z)<1e-10);
val = besselj(n,z)./(z).^n;
%val = mfun('BesselJ',n,z)./(z).^n; %->SLOW
val(ind_small) = 2^(-n)/factorial(n)*(1-(z(ind_small)/2).^2/(n+1));
warning('ON','all');
end