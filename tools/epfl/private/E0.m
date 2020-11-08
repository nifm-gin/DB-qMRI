%% [m,param] = E0(x, param)
%
% Function that performs the matrix-vector multiplication that models a 2D MRI
% experiment with a uniform receiving coil sensitivity and a given arbitrary
% k-space trajectory.
%
% INPUTS:       * the original 2D image
%               * a structure containing the parameters (k-space trajectory,
% 			method for computation, ...).
%
% OUTPUTS:      * the k-space data vector
%               * (optionally) a structure carrying precomputed data to fasten
% 			future calls.
%
% SEE: EH0.m EHE0.m
%
% Copyright, Matthieu Guerquin-Kern, 2012

function [m,param] = E0(x, param)

res = size(x);
param.res = res;

if ~isfield(param,'method')
    if max(res)>32
        param.method='nuft gg';
    else
        param.method='dft mex';
    end
end

switch lower(param.method)
    case 'matlab script'
        % Good as a reference for the other routines. The complexity is N^4.
        % It is horribly slow for 2D images of dimension larger than N=100.
        Nsamples = size(param.k,1);
        m = zeros(Nsamples,1);
        [X,Y] = ndgrid((0:res(1)-1)/res(1),(0:res(2)-1)/res(2));
        for i=1:Nsamples
            m(i) = x(:).'*exp(1i*2*pi*(param.k(i,1)*X(:)+param.k(i,2)*Y(:)));
        end
    case 'dft mex'
        % efficient MEX implementation of 'matlab script'. Seems to work
        % fine and shows very good accuracy in our tests. The code is
        % optimized, hence complicated, so it might be buggy.
        % The complexity is N^4. For large dataset, prefer the two next options.
        m = nuft_forw(x,param.k(:,1),param.k(:,2));
    case 'nuft gg'
        % MEX implementation of Greengard's fast gaussian gridding.
        % Complexity N^2log(N).
        if ~isfield(param,'st')
            param.st = nuft_gg_init(-2*pi*[param.k(:,1)/res(1),param.k(:,2)/res(2)], res, 12, 4*res);
        end
        m = nuft_gg_forw(x, param.st);
    case 'nuft gg fast'
        % MEX implementation of Greengard's fast gaussian gridding.
        % Complexity N^2log(N).
        if ~isfield(param,'st')
            param.st = nuft_gg_init(-2*pi*[param.k(:,1)/res(1),param.k(:,2)/res(2)], res, 8, 3*res);
        end
        m = nuft_gg_forw(x, param.st);
    case 'nufft'
        % NUFFT implementation by Fessler etal. Poor accuracy, minimum
        % 1e-10, particularly bad for small dimensions. Very fast at large
        % scale. For very large scale, it causes memory problems.
        % The corresponding code is not packaged here. You will find it at
        % http://web.eecs.umich.edu/~fessler/irt/irt/nufft .
        if ~exist('nufft','file')
            error('The IRT package seems not installed.\nYou will find it at:\n http://web.eecs.umich.edu/~fessler/irt/irt/nufft');
        end
        if ~isfield(param,'st')
            param.st = nufft_init(-2*pi*[param.k(:,1)/res(1),param.k(:,2)/res(2)], res, 10*[1,1], 2*res);
        end
        m = nufft(x, param.st);
    otherwise
        error(sprintf('unknown method ''%s''',param.method));
end

if isfield(param,'weight')
    m = m.*param.weight;
end
