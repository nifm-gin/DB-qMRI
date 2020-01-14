function [x,x_iter] = conjugate_gradient(E,b,tol,maxit,x,verbose,newfigure)
% A conjugate gradient solver that allows to plot the result in each
% iteration
% Solves: E*x=b
%
% x = conjugate_gradient(E,b)
% x = conjugate_gradient(E,b,tol)
% x = conjugate_gradient(E,b,tol,maxit)
% x = conjugate_gradient(E,b,tol,maxit,z)
% x = conjugate_gradient(E,b,tol,maxit,z,verbose)
% x = conjugate_gradient(E,b,tol,maxit,z,verbose,newfigure)
% [x,x_iter] = conjugate_gradient(_____)
%
%
% Input:
%   E         =  Hermitian and positive definite matrix or operator that
%                immitates such a matrix and implements a mtimes funciton.
%                Alternatively, A can function handle that imitates such a
%                matrix.
%   b         =  Right hand side of the equation A*x=b
%   tol       =  Tolerance of the method (default: 1e-6)
%   maxit     =  maximum number of iterations (default: 20)
%   x         =  initial guess (default: 0)
%   verbose   =  0 for no output, 1 for plotting the images and in each
%                iteration and 2 for printing the residal of the for each
%                iteration.
%   newfigure =  0 for using the same figure for each call of the funciton
%                or 1 for opening a new figure.
%
%
% Output:
%   x      = Result
%   x_iter = Result after each iteration
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) Jakob Asslaender, August 2016
% New York University School of Medicine, Center for Biomedical Imaging
% University Medical Center Freiburg, Medical Physics
% jakob.asslaender@nyumc.org
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


persistent h;
if isa(E, 'function_handle')
    function_flag = 1;
else
    function_flag = 0;
end

if nargin<=2 || isempty(tol)
    tol = 1e-6;
end
if nargin<=3 || isempty(maxit)
    maxit = 20;
end
if nargin<=4 || isempty(x)
    if isnumeric(b)
        x = zeros(size(b));
    else
        eval(['z = ',class(b),'(size(b));']);
    end
    d = b;
elseif isnumeric(x) && isscalar(x) && x==0
    x = zeros(size(b));
    d = b;
else
    if function_flag
        d = b - E(x);
    else
        d = b - E*x;
    end
end
if nargin<=5
    verbose = 1;
end
if nargin<=6
    newfigure=0;
end

r = d;

normrr0 = sum(col(conj(b) .* b));
normrr  = sum(col(conj(r) .* r));

for n=1:maxit
    if function_flag
        Ed = E(d);
    else
        Ed = E*d;
    end        
    tmp = conj(d) .* Ed;
    alpha = normrr/sum(col(tmp));
    x = x + alpha*d;
    r = r - alpha*Ed;
    normrr2 = sum(col(conj(r) .* r));
	beta = normrr2/normrr;
    normrr = normrr2;
	d = r + beta*d;
    
    if nargout==2
        x_iter{n} = x;
    end
    
    if verbose==1
        if isempty(h) || ~ishandle(h) || ~strcmp(get(h,'Tag'),'cg_figure') || (newfigure && n==1);
            h = figure('Tag','cg_figure');
        end
        if gcf~=h
            set(0,'CurrentFigure',h);
        end
        if isnumeric(x) && isvector(x)
            plot(abs(x));
        else
            imagesc34d(abs(x),1); colormap gray;
        end
        axis off;
        title(['CG iteration #', num2str(n),': r/r0 = ', num2str(sqrt(normrr/normrr0))]);
        drawnow;
    elseif verbose==2
        fprintf(['CG iteration #', num2str(n),': r/r0 = ', num2str(sqrt(normrr/normrr0)), '\n']);
    end
    
    if sqrt(normrr/normrr0) < tol
        break;
    end
    
end
