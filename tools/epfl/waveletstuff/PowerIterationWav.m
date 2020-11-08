% PowerIterationWav.m
%
% Computes the greatest eigenvalue alpha of a square matrix M.
%            for all x,  ||Mx||_2 <= alpha*||x||_2
%
% If M is positive definite and expressed as M = AT*A, we have
%            for all x,  ||Ax||^2 <= alpha*||x||^2,
% and alpha can be seen as the Lipschitz constant of the function f(x)= Ax.
%
% INPUT:        * M as a matrix or function handle
%               * v an initial random vector
%
% OUTPUT:       * alpha the greatest (in absolute value) eigenvalue
%               * v the associated eigenvector
%
% Matthieu Guerquin-Kern, Biomedical Imaging Group / EPF Lausanne,
% 11-08-2009 (dd-mm-yyyy)

function alpha = PowerIterationWav(M,w)

if isa(w,'double')
    [alpha,v] = power_iteration(M,w);
elseif ~isa(w,'wavelet')
    error('unknown type of entry!');
end

if ~isa(M, 'function_handle')
    if nargin<2
        w = rand(size(M,2),1);
    end
    M = @(x) M*x;
end

tol = 1e-3;
maxit = 5e1;
J = sort(w.depth);
Jm = circshift(J(:),-1);Jm(1) = 0;
m = 1 + sum(double(J(:)-Jm(:)).*(2.^(numel(J):-1:1)'-1));clear Jm;
w0 = 0*w;
alpha = zeros(m,m);
for i = 1:m       
    for j = 1:i
        %disp(sprintf('i: %d; j: %d\n',i,j));
        tmp = w{i};wi = w0;
        wi{i} = tmp;
        wi = wi/norm(wi(:),2);
        Mwi = M_ij(wi,M,j,i);
        crit = 1/eps;
        counter=0;
        while counter<maxit && crit>tol
            tmp = Mwi(:)'*wi(:);
            if norm(real(tmp(:)))<1e5*norm(imag(tmp(:)))
                warning('power_iteration_wavelets: the matrix provided is not positive definite!');
            end
            alpha(i,j) = real(tmp);
            %sum(sqrt(alpha),1)
            wi = Mwi/norm(Mwi(:),2);
            Mwi = M_ij(wi,M,j,i);
            diff = Mwi-alpha(i,j)*wi;
            crit = norm(diff(:),2);
            counter = counter+1;
        end
        if i~=j, alpha(j,i) = alpha(i,j);end
        %disp(sprintf('counter: %d; criterion: %g\n',counter,crit));
    end
end
alpha = sqrt(abs(alpha));
alpha = sum(alpha,2);
end

function wi = M_ij(wi,M,j,i)
    wj = M(wi);   tmp = wj{j};
    wj = 0*wj;    wj{j} = tmp;
    wi = M(wj);   tmp = wi{i};
    wi = 0*wi;    wi{i} = tmp;
end