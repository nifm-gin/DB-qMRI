%% MRDataBezier.m
%
% Returns the MR data corresponding to the indicator function of a region
% whose contour is defined by a Bezier curve of order 1 (i.e. a polygon) or
% 2, that is modulated by a polynomial or sinusoidal profile, for the
% given kspace samples.
%
% Matthieu Guerquin-Kern, Biomedical Imaging Group / EPF Lausanne,
% 15-02-2011 (dd-mm-yyyy)

function m = MRDataBezier(region,s,k)

weight = region.weight;
switch lower(region.type)
    case 'polygon'
        r = region.vertex;
        beta = circshift(r,[-1,0])-r;
        gamma = zeros(size(beta));
    case 'bezier'
        r = control2node(region.control);
        c = circshift(region.control,[1,0]);
        rp1 = circshift(r,[-1,0]);
        beta = c-circshift(c,[1,0]);
        gamma = rp1+r-2*c;
end

switch lower(s.model)
    case 'polynomial'
        m = weight*MRDataBezierPolynomial(s.data,k,r,beta,gamma);
    case 'sinusoidal'
        m = weight*MRDataBezierSinusoidal(s.data,k,r,beta,gamma);
    otherwise
        error('unknown sensitivity model. Must be either ''polynomial'' or ''sinusoidal''.')
end
end

%%
function m = MRDataBezierSinusoidal(s,k,r,beta,gamma)
Ns = length(s);
%Nk = size(k,2);
L = floor(sqrt(Ns));
[kx,ky] = GenerateFullCart2DGrid(L*[1,1]);
kx = kx(:)/2;
ky = ky(:)/2;
Nmax = 1e7; % ADJUST THIS VALUE DEPENDING ON THE AMOUNT OF RAM AVAILABLE (IN ORDER TO AVOID SWAPING)
Nnodes = size(beta,1);
Nsamples = size(k,2);
m = zeros(1,Nsamples);
N = Nsamples*Nnodes;
n = floor(N/Nmax);
n0 = floor(Nsamples/max(1,n));
for c=1:length(kx)
    kn = [k(1,:)+kx(c);k(2,:)+ky(c)];
    for i = 0:n-1
        f = calculatef(r,beta,gamma,-2*pi*kn(:,i*n0+1:(i+1)*n0),0);
        m(i*n0+1:(i+1)*n0) = m(i*n0+1:(i+1)*n0) + s(c)*f{1,1};
    end
    if Nsamples>n*n0
        f = calculatef(r,beta,gamma,-2*pi*kn(:,n*n0+1:end),0);
        m(n*n0+1:end) = m(n*n0+1:end) + s(c)*f{1,1};
    end
end
end

function m = MRDataBezierPolynomial(s,k,r,beta,gamma)
% Laurent Leujeune, semester project, BIG EPFL, winter 2010-2011
% modified by MGK, march 2011
N = length(s);
D =  floor((sqrt(1+N*8)-3)/2); % works only in 2D!
f = calculatef(r,beta,gamma,-2*pi*k,D);
[a1 a2] = meshgrid(0:D, 0:D);
m = zeros(1,size(k,2));
for d=0:D
    ind = find(a1 + a2 <= d);
    alpha = [a1(ind) a2(ind)]';
    [tmp, I] = sort(sum(alpha));
    alpha = [alpha(1,I); alpha(2,I)];
    for j=(d+1)*(d+2)/2-d:(d+1)*(d+2)/2
        m = m + s(j)*f{alpha(1,j)+1, alpha(2,j)+1};
    end
end
end

function f = calculatef(r,beta,gamma,w,D)
%UNTITLED Summary of this function goes here
%   f: 1st dimension vector alpha, 2nd and 3rd: components of
% fprintf('calculatef...\n');
%
% Implements Theorems 2.4 and 2.5.
%
% Laurent Leujeune, semester project, BIG EPFL, winter 2010-2011
% modified by MGK, march 2011

[a1 a2] = meshgrid(0:D, 0:D);
ind = find(a1 + a2 <= D); % find the set of vectors alpha s.t. |alpha|<=D
alpha = [a1(ind) a2(ind)]';
[tmp, I] = sort(sum(alpha));
alpha = flipud([alpha(1,I); alpha(2,I)]);

Nw = size(w,2);
normwsquared = sum(w.^2,1);
indknull = find(sqrt(normwsquared)<=10*eps); % locate the case w=0 for separate treatment
indknotnull = find(sqrt(normwsquared)>10*eps);
normwsquared = normwsquared(indknotnull);
w = w(:,indknotnull);
g = calculateg(r,beta,gamma,w,alpha);
g0 = calculateg0(r,beta,gamma, alpha);

N = size(r,1); % number of control points or vertex
f = cell(D+1, D+1);

for i=1:size(alpha,2) % for each alpha s.t. |alpha|<=D
    %fprintf('alpha= (%d %d)\n',alpha(1,i),alpha(2,i));
    curr_alpha = repmat(alpha(:,i),1,N);
    p = [repmat(0:curr_alpha(1,1),1,curr_alpha(2,1)+1);sort(repmat(0:curr_alpha(2,1),1,curr_alpha(1,1)+1))];
    calp1 = curr_alpha(1,1)+1;
    calp2 = curr_alpha(2,1)+1;
    f{calp1,calp2} = zeros(1,Nw);
    for n=1:size(p,2)
        curr_p = repmat(p(:,n),1,size(w,2));
        curr_alpha_minus_p = repmat(curr_alpha(:,1)-curr_p(:,1),1,size(w,2));
        f{calp1,calp2}(indknotnull) = f{calp1,calp2}(indknotnull) + prod((-1i*w./repmat(normwsquared,2,1)).^curr_alpha_minus_p,1)...
            *factorial(sum(curr_alpha_minus_p(:,1),1))*...
            nchoosekMultiInd(curr_alpha(:,1),curr_p(:,1)).*transpose(g{curr_p(1,1)+1,curr_p(2,1)+1});
    end
    f{calp1,calp2}(indknotnull) = (1i./normwsquared).*f{calp1,calp2}(indknotnull);
    if ~isempty(indknull)
        f{calp1,calp2}(indknull) = g0(calp1,calp2);
    end
end
end

function g = calculateg(r,beta,gamma, w, alpha)
%Builds matrices g for each value of alpha
%INPUTS:
%       element: - Polygon case
%                - ... or Bezier curve structure
%       k: k-space trajectory (2xN)
%       alpha: each row represents a vector alpha
% Laurent Leujeune, semester project, BIG EPFL, winter 2010-2011
% modified by MGK, march 2011

N = size(r,1);
g = cell(max(alpha(1,:))+1, max(alpha(2,:))+1 );
leftterm = exp(-1i*w'*r');
beta_normal = beta*[0 -1; 1 0];
gamma_normal = gamma*[0 -1; 1 0];
orders = 0:(1+2*max(sum(alpha,1)));

h = buildH(orders,beta*w,gamma*w);
for i=1:size(alpha,2)
    curr_alpha = alpha(:,i);
    calp1 = curr_alpha(1,1)+1;
    calp2 = curr_alpha(2,1)+1;
    [L,indL] = findCombinationsForL( curr_alpha );
    rightterm = zeros(N,size(w,2));
    for n=1:size(indL,2)
        l1 = repmat(L(:,indL(1,n))',N,1);
        l2 = repmat(L(:,indL(2,n))',N,1);
        l3 = repmat(L(:,indL(3,n))',N,1);
        rightterm = rightterm...
            + prod(factorial(curr_alpha))/(prod(factorial(l1(1,:)))*prod(factorial(l2(1,:)))*prod(factorial(l3(1,:))))...
            *repmat( prod(r.^l1,2).*prod(beta.^l2,2).*prod(gamma.^l3,2) ,1,size(w,2)).*...
            (  (beta_normal*w).*h{1 + sum(l2(1,:)) + 2*sum(l3(1,:))} + ...
            2*(gamma_normal*w).*h{2 + sum(l2(1,:)) + 2*sum(l3(1,:))} );
    end
    clear curr_alpha;
    g{calp1,calp2} = sum(leftterm.*rightterm.',2);
end
end

function g = calculateg0(r,beta,gamma, alpha)
% This function computes the alpha-th order moment of the region.
% For alpha=[0,0] it corresponds to the area of the region.
%
% Laurent Leujeune, semester project, BIG EPFL, winter 2010-2011
% modified by MGK, march 2011

N = size(r,1);
ek = [1;0];
g = zeros(max(alpha(1,:))+1, max(alpha(2,:))+1 );
beta_normal = beta*[0 -1; 1 0];
gamma_normal = gamma*[0 -1; 1 0];

orders = 0:(1+2*(1+max(sum(alpha,1))));
h = 1./(orders+1);
alpha_plus_ek = alpha + repmat(ek,1,size(alpha,2));

for i=1:size(alpha,2)
    curr_alpha = alpha(:,i);
    curr_alpha_plus_ek = repmat(alpha_plus_ek(:,i),1,N);
    tmp = 1/curr_alpha_plus_ek(1,1);
    calp1 = curr_alpha(1,1)+1;
    calp2 = curr_alpha(2,1)+1;
    [L,indL] = findCombinationsForL( curr_alpha_plus_ek(:,1) );
    rightterm = zeros(N,1);
    for k=1:size(indL,2)
        l1 = repmat(L(:,indL(1,k))',N,1);
        l2 = repmat(L(:,indL(2,k))',N,1);
        l3 = repmat(L(:,indL(3,k))',N,1);
        rightterm = rightterm...
            + (tmp)*prod(r.^l1,2).*prod(beta.^l2,2).*prod(gamma.^l3,2).*...
            prod(factorial(curr_alpha_plus_ek(:,1)))/(prod(factorial(l1(:,1)))*prod(factorial(l2(:,1)))*prod(factorial(l3(:,1)))).*...
            (  (beta_normal*ek).*h(1 + sum(l2(1,:)) + 2*sum(l3(1,:))) + ...
            2*(gamma_normal*ek).*h(2 + sum(l2(1,:)) + 2*sum(l3(1,:)))      );
    end
    g(calp1,calp2) = sum(rightterm);
end
end

function h = buildH( orders, a, b)
%	hlayers: cells where the first cell is the function h for m = 0
%	orders: vector containing the orders at which the derivatives have to be calculated: [0 ... |p|max].
%	w: k-space coordinates
% Laurent Leujeune, semester project, BIG EPFL, winter 2010-2011
% modified by MGK, march 2011

sz = size(a);

%fprintf('\nBUILDH for %d elements\n\n',numel(a));
m_max = max(orders);
h = cell(1,m_max+1);
for m = 0:m_max
    h{m+1} = zeros(sz);
end

%fprintf('statistics: min a: %f, max a: %f, min b: %f, max b: %f\n',min(a(:)),max(a(:)),min(b(:)),max(b(:)));
%fprintf('statistics: std a: %f, std b: %f\n',std(a(:)),std(b(:)));
%figure;plot(a(:),b(:),'.');drawnow

za_ind = (abs(a) <= 2); %index of small a
zb_ind = (abs(b) <= 1e-5); %index of small b
zazb_ind = find(za_ind&zb_ind); % -> use backwards
nzazb_ind = find((~za_ind)&zb_ind); % -> use series
nzb_ind = find(~zb_ind); % -> use forward
clear za_ind zb_ind;

b_zazb = b(zazb_ind);
a_zazb = a(zazb_ind);
b_nzazb = b(nzazb_ind);
a_nzazb = a(nzazb_ind);
b_nzb = b(nzb_ind);
a_nzb = a(nzb_ind);
clear a b;

%*************************************************************************
% Lets deal first with the case a and b small
% the backward iteration is stable in that case
if numel(a_zazb)>0
    %fprintf('Using backward recursion. Target order: %d, Max a %g, Min a %g, Max b %g, Min b %g\n',m_max,max(abs(a_zazb(:))),min(abs(a_zazb(:))),max(abs(b_zazb(:))),min(abs(b_zazb(:))));
    h = BuildH_backward(a_zazb,b_zazb,h,zazb_ind);
end

%*************************************************************************
% Second we deal with a large and b small
% in that case the series development in b are stable
if numel(a_nzazb)>0
    %fprintf('Using series. Target order %d\n',m_max);
    h = BuildH_series(a_nzazb,b_nzazb,h,nzazb_ind);
end

%*************************************************************************
% Third, we consider the case where b is large
% in that case, the forward iteration is adequate
if numel(a_nzb)>0
    %fprintf('Using forward recursion. Target order %d\n',m_max);
    h = BuildH_forward(a_nzb,b_nzb,h,nzb_ind);
end
end

function h = BuildH_forward(a,b,h,index)
m_max = length(h)-1;
implementation = 'fastest';
e =     myerfz(sqrt(2)./sqrt(abs(b))/4.*real(a+2*b),implementation)...
    -  myerfz(sqrt(2)./sqrt(abs(b))/4.*real(a),implementation);
e(b<0) = conj(e(b<0));
h{1}(index) = sqrt(pi)*exp(1i*a.^2./(4*b))./(2*sqrt(1i)*sqrt(b)).*e;clear e;
diverged = (abs(h{1}(index))>1);
if any(diverged(:))
    %fprintf('forward iteration diverged at order 0\n');
    figure(99);hold on;plot(a(diverged),b(diverged),'.r');hold off;xlabel('a');ylabel('b');title(sprintf('places where the computation of h^{(%d)} diverged',0));
end
if m_max>=1
    Eab = exp(-1i*(a+b));
    h{2}(index) = (-a.*h{1}(index) + 1j*Eab - 1j)./(2*b);
    diverged = (abs(h{2}(index))>0.5);
    if any(diverged(:))
        %fprintf('forward iteration diverged at order %d\n',1);
        figure(99);hold on;plot(a(diverged),b(diverged),'.r');hold off;xlabel('a');ylabel('b');title(sprintf('places where the computation of h^{(%d)} diverged',1));
    end
    for mk=3:m_max+1
        h{mk}(index) = (-1i*a.*h{mk-1}(index)  + (mk-2)*h{mk-2}(index) - Eab)./(1i*2*b);
        diverged = (mk*abs(h{mk}(index))>1);
        if any(diverged(:))
            %fprintf('forward iteration diverged at order %d\n',mk-1);
            figure(99);hold on;plot(a(diverged),b(diverged),'.r');hold off;xlabel('a');ylabel('b');title(sprintf('places where the computation of h^{(%d)} diverged',mk-1));
        end
    end
end
end

function h = BuildH_series(a,b,h,index)
m_max = length(h)-1;
% We build the table of h^(mk)(a,0)
Nseries=1;
ha0 = cell(1,2*Nseries+m_max);
ha0{1} = exp(-1i*a/2).*sinc(a/(2*pi));
Ea = exp(-1i*a);
delta = eps*ones(size(Ea));
for m = 1:(2*Nseries+m_max)
    ha0{m+1} = (m*ha0{m}-Ea)./(1j*a);
    delta = m*delta./abs(a);
    diverged = ((m+1)*abs(ha0{m+1})>1);
    if any(diverged(:))
        %fprintf('series: ha0: divergence at order %d\n',m);
    end
end
%fprintf('\t-> a large and b small (table). Accuracy: %g\n',max(delta(:)./abs(ha0{m+1})));
% We build the series out of it, for each m between 0 and m_max
for m = 0:m_max
    %h{m+1}(index) = 0;
    c = ones(size(b));
    for n=0:Nseries
        h{m+1}(index) = h{m+1}(index) + c.*ha0{1+m+2*n};%(-1j*b).^i.*ha0{1+m+2*i}/factorial(i);%
        c = c.*(-1j*b)/(n+1);
    end
    diverged = ((m+1)*abs(h{m+1}(index))>1);
    if any(diverged(:))
        %fprintf('series: divergence at order %d\n',m);
        figure(99);hold on;plot(a(diverged),b(diverged),'.g');hold off;xlabel('a');ylabel('b');title(sprintf('places where the computation of h^{(%d)}(a,b) diverged',m));
    end
end
end

function h = BuildH_backward(a,b,h,index)
m_max = length(h)-1;
m_start = m_max+40;
% initialize
hb = cell(1,m_max+1);
for m = 1:m_start
    hb{m+1} = zeros(size(a));
end

Eab = exp(-1i*(a+b));
hb{m_start+1} = Eab/(m_start+1);
hb{m_start} = ( 1j*a.*hb{m_start+1} + Eab)/m_start;
delta = ones(size(Eab))/(m_start+1);
delta_old = ones(size(Eab))/(m_start+2);
for mk=m_start-1:-1:m_max+2
    hb{mk} = ( 1i*a.*hb{mk+1} + 2*1i*b.*hb{mk+2} + Eab)/mk;
    dt = abs(a.*delta+2*b.*delta_old)/mk;
    delta_old = delta;
    delta = dt;
end
for mk=m_max+1:-1:1
    hb{mk} = ( 1i*a.*hb{mk+1} + 2*1i*b.*hb{mk+2} + Eab)/mk;
    dt = abs(a.*delta+2*b.*delta_old)/mk;
    delta_old = delta;
    delta = dt;
    diverged = (mk*abs(hb{mk})>1);
    if any(diverged(:))
        %fprintf('backward: divergence at order %d\n',mk-1);
        figure(99);hold on;plot(a(diverged),b(diverged),'.b');hold off;xlabel('a');ylabel('b');title(sprintf('places where the computation of h^{(%d)}(a,b) diverged',mk-1));
    end
end
for i=1:m_max+1
    h{i}(index) = hb{i};
end
% controlling the result at order 0
e =     myerfz(sqrt(2)./sqrt(abs(b))/4.*real(a+2*b),'fastest')...
    -  myerfz(sqrt(2)./sqrt(abs(b))/4.*real(a),'fastest');
e(b<0) = conj(e(b<0));
h_control = sqrt(pi)*exp(1i*a.^2./(4*b))./(2*sqrt(1i)*sqrt(b)).*e;clear e;
h_control(abs(b)<50*eps) = exp(-1i*a(abs(b)<50*eps)/2).*sinc(a(abs(b)<50*eps)/(2*pi));
%fprintf('\t-> a and b small. Accuracy: %g and %g\n',max(delta(:)./abs(hb{1})),max(abs(hb{1}-h_control(:))./abs(hb{1})));
end

function [L,indL] = findCombinationsForL( alpha )
%Used to compute g (bezier case). We need to find linear combinations
% of 3 vectors l1,l2,l3 where li can be any vector in L= [0 ... alpha]
%and must respect the constraint: l1 + l2 + l3 = alpha
%
% OUTPUT:   L= [0 ... alpha]
%           indL: 3x(2*|alpha| + 1)
%                   one column per combination
%                   rows indicate index of vector from L
%
% Laurent Leujeune, semester project, BIG EPFL, winter 2010-2011
D = sum(alpha);
[a1 a2] = meshgrid(0:D, 0:D);
ind = find(a1 + a2 <= D);
L = [a1(ind) a2(ind)]';

ind = npermutek(1:size(L,2),3)';

R = L(:,ind(1,:)) + L(:,ind(2,:)) + L(:,ind(3,:));
alpha = repmat(alpha,1,size(R,2));

indL = ind(:,(0.5*sum(R == alpha) == 1));

end

function [Matrix,Index] = npermutek(N,K)
%NPERMUTEK Permutation of elements with replacement/repetition.
% Author:  Matt Fig
% Contact:  popkenai@yahoo.com
% Copyright (c) 2006, Matt Fig
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without 
% modification, are permitted provided that the following conditions are 
% met:
% 
%     * Redistributions of source code must retain the above copyright 
%       notice, this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright 
%       notice, this list of conditions and the following disclaimer in 
%       the documentation and/or other materials provided with the distribution
%       
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE 
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF 
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS 
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN 
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) 
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE 
% POSSIBILITY OF SUCH DAMAGE.

if nargin ~= 2
    error('NPERMUTEK requires two arguments. See help.')
end

if isempty(N) || K == 0,
    Matrix = [];
    Index = Matrix;
    return
elseif floor(K) ~= K || K<0 || ~isreal(K) || numel(K)~=1
    error('Second argument should be a real positive integer. See help.')
end

LN = numel(N);  % Used in calculating the Matrix and Index.

if K==1
    Matrix = N(:); % This one is easy to calculate.
    Index = (1:LN).';
    return
elseif LN==1
    Index = ones(K,1);
    Matrix = N(1,Index);
    return
end

CLS = class(N);

if ischar(N)
    CLS = 'double';  % We will deal with this at the end.
end

L = LN^K;  % This is the number of rows the outputs will have.
Matrix = zeros(L,K,CLS);  % Preallocation.
D = diff(N(1:LN));  % Use this for cumsumming later.
LD = length(D);  % See comment on LN.
VL = [-sum(D) D].';  % These values will be put into Matrix.
% Now start building the matrix.
TMP = VL(:,ones(L/LN,1,CLS));  % Instead of repmatting.
Matrix(:,K) = TMP(:);  % We don't need to do two these in loop.
Matrix(1:LN^(K-1):L,1) = VL;  % The first column is the simplest.

if nargout==1
    % Here we only have to build Matrix the rest of the way.
    for ii = 2:K-1
        ROWS = 1:LN^(ii-1):L;  % Indices into the rows for this col.
        TMP = VL(:,ones(length(ROWS)/(LD+1),1,CLS));  % Match dimension.
        Matrix(ROWS,K-ii+1) = TMP(:);  % Build it up, insert values.
    end
    
else
    % Here we have to finish Matrix and build Index.
    Index = zeros(L,K,CLS);  % Preallocation.
    VL2 = ones(size(VL),CLS);  % Follow the logic in VL above.
    VL2(1) = 1-LN;  % These are the drops for cumsum.
    TMP2 = VL2(:,ones(L/LN,1,CLS));  % Instead of repmatting.
    Index(:,K) = TMP2(:);  % We don't need to do two these in loop.
    Index(1:LN^(K-1):L,1) = 1;
    
    for ii = 2:K-1
        ROWS = 1:LN^(ii-1):L;  % Indices into the rows for this col.
        F = ones(length(ROWS)/(LD+1),1,CLS);  % Don't do it twice!
        TMP = VL(:,F);  % Match dimensions.
        TMP2 = VL2(:,F);
        Matrix(ROWS,K-ii+1) = TMP(:); % Build them up, insert values.
        Index(ROWS,K-ii+1) = TMP2(:);
    end
    
    Index(1,:) = 1;  % The first row must be 1 for proper cumsumming.
    Index = cumsum(Index);  % This is the time hog.
end

Matrix(1,:) = N(1);  % For proper cumsumming.
Matrix = cumsum(Matrix);  % This is the time hog.

if ischar(N)
    Matrix = char(Matrix);  % char was implicitly cast to double above.
end

end

function [ out ] = nchoosekMultiInd( p, q )
%   Product of the number of combinations of pi things taken qi at a time.
% Laurent Leujeune, semester project, BIG EPFL, winter 2010-2011
% modified by MGK, april 2011

if any(min(p-q)<0)
    error('myApp:argChk', 'components pi must be greater than qi');
end

out = prod(factorial(p))./( prod(factorial(p-q)).*prod(factorial(q)));
end
