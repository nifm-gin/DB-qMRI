%% [x,param] = EH0(m, param)
%
% Function that performs the matrix-vector multiplication with the Hermitian
% transpose of the matrix that models a 2D MRI experiment with a uniform
% receiving coil sensitivity and a given arbitrary k-space trajectory.
%
% INPUTS:       * the k-space data vector
%               * a structure containing the parameters (k-space trajectory,
% 			method for computation, ...).
%
% OUTPUTS:      * the backprojected 2D image
%               * (optionally) a structure carrying precomputed data to fasten
% 			future calls.
%
% SEE: E0.m EHE0.m
%
% Copyright, Matthieu Guerquin-Kern, 2012

function [x,param] = EH0(m, param)

if ~isfield(param,'res')
    error('the image dimensions must be precised for computation');
else
    res = param.res;
end
if isfield(param,'weight')
    m = m.*conj(param.weight);
end
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
        x = zeros(res);
        [X,Y] = ndgrid((0:res(1)-1)/res(1),(0:res(2)-1)/res(2));
        for i=1:Nsamples
            x = x+m(i).'*exp(-1i*2*pi*(param.k(i,1)*X+param.k(i,2)*Y));
        end
    case 'dft mex'
        % efficient MEX implementation of 'matlab script'. Seems to work
        % fine and shows very good accuracy in our tests. The code is
        % optimized, hence complicated, so it might be buggy.
        % The complexity is N^4. For large dataset, prefer the two next options.
        x = nuft_back(m,param.k(:,1),param.k(:,2),res);
    case 'nuft gg'
        % MEX implementation of Greengard's fast gaussian gridding.
        % Complexity N^2log(N).
        if ~isfield(param,'st')
            param.st = nuft_gg_init(-2*pi*[param.k(:,1)/res(1),param.k(:,2)/res(2)],res,12,4*res);
        end
        x = nuft_gg_back(m, param.st);
    case 'nuft gg fast'
        % MEX implementation of Greengard's fast gaussian gridding.
        % Complexity N^2log(N).
        if ~isfield(param,'st')
            param.st = nuft_gg_init(-2*pi*[param.k(:,1)/res(1),param.k(:,2)/res(2)],res,8,3*res);
        end
        x = nuft_gg_back(m, param.st);
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
            st = nufft_init(-2*pi*[param.k(:,1)/res(1),param.k(:,2)/res(2)], res, [10,10], 2*res);
        end
        x = nufft_adj(m, param.st);
    otherwise
	error(sprintf('unknown method ''%s''',param.method));
end

x = x/prod(res);
