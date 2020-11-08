% PowerIteration.m
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

function [alpha,v] = PowerIteration(M,v)

if ~isa(M, 'function_handle')
    if nargin<2
        v = rand(size(M,2),1);
    end
    M = @(x) M*x;
end
v = v/norm(v(:),2);
Mv = M(v);
alpha = Mv(:)'*v(:);
diff = Mv-alpha*v;
crit = norm(diff(:),2);
tol = 1e-3;
maxit = 5e1;
counter=0;

while counter<maxit && crit>tol
    alpha = Mv(:)'*v(:);
    v = Mv/norm(Mv(:),2);
    Mv = M(v);
    diff = Mv-alpha*v;
    crit = norm(diff(:),2);
    counter = counter+1;
end
alpha = real(alpha);