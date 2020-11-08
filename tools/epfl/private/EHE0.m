%% [y,param] = EHE0(x, param)
%
% Function that performs the matrix-vector multiplication E0H*E0 with E0 the
% matrix that models a 2D MRI experiment with a uniform receiving coil
% sensitivity and a given arbitrary k-space trajectory and E0H it Hermitian
% transpose.
% This implementation is equivalent but faster than the consecutive use of
% functions E0 and EH0.
%
% INPUTS:       * the original 2D image
%               * a structure containing the parameters (kernel,...).
%
% OUTPUTS:      * the resulting 2D image
%               * (optionally) a structure carrying precomputed data to fasten 
%			future calls.
%
% SEE: E0.m EH0.m
%
% Copyright, Matthieu Guerquin-Kern, 2012

function [y,param] = EHE0(x, param)

res = size(x);

if ~isfield(param,'kernel')
    if ~isfield(param,'method')
        if max(res)>32
            param.method='nuft gg';
        else
            param.method='dft mex';
        end
    end
    if isfield(param,'weight')
        w=abs(param.weight).^2;
    else
        w=ones(size(param.k,1),1);
    end
    t0 = clock();
    param.kernel = zeros(2*res);
    switch lower(param.method)
        case 'matlab script'
            % Good as a reference for the other routines. The complexity is N^4.
            % It is horribly slow for 2D images of dimension larger than N=100.
            kernel = zeros(2*res-1);
            [X,Y] = ndgrid((-res(1)+1:res(1)-1)/res(1),(-res(2)+1:res(2)-1)/res(2));
            for i=1:size(param.k,1)
                kernel = kernel+w(i)*exp(-1i*2*pi*(param.k(i,1)*X+param.k(i,2)*Y));
            end
            param.kernel(2:end,2:end) = kernel;
        case 'dft mex'
            % efficient MEX implementation of 'matlab script'. Seems to work
            % fine and shows very good accuracy in our tests. The code is
            % optimized, hence complicated, so it might be buggy.
            % The complexity is N^4. For large dataset, prefer the two next options.
            kernel = nuft_kern(param.k(:,1),param.k(:,2), res, w);
            param.kernel(2:end,2:end) = kernel;
        case 'nuft gg'
            % MEX implementation of Greengard's fast gaussian gridding.
            % Complexity N^2log(N).
            st = nuft_gg_init(-2*pi*[param.k(:,1)/res(1),param.k(:,2)/res(2)],2*res,12,8*res);
            param.kernel = nuft_gg_back(w.*exp(2*pi*1j*(param.k(:,1)+param.k(:,2))), st);
        case 'nuft gg fast'
            % MEX implementation of Greengard's fast gaussian gridding.
            % Complexity N^2log(N).
            st = nuft_gg_init(-2*pi*[param.k(:,1)/res(1),param.k(:,2)/res(2)],2*res,8,6*res);
            param.kernel = nuft_gg_back(w.*exp(2*pi*1j*(param.k(:,1)+param.k(:,2))), st);
        case 'nufft'
            % NUFFT implementation by Fessler etal. Poor accuracy, minimum
            % 1e-10, particularly bad for small dimensions. Very fast at large
            % scale. For very large scale, it causes memory problems.
            % The corresponding code is not packaged here. You will find it at
            % http://web.eecs.umich.edu/~fessler/irt/irt/nufft .
            if ~exist('nufft','file')
                error('The IRT package seems not installed.\nYou will find it at:\n http://web.eecs.umich.edu/~fessler/irt/irt/nufft');
            end
            st = nufft_init(-2*pi*[param.k(:,1)/res(1),param.k(:,2)/res(2)], 2*res, 10*[1,1], 4*res);
            param.kernel = nufft_adj(w.*exp(2*pi*1j*(param.k(:,1)+param.k(:,2))),st);
        otherwise
            error(sprintf('unknown method ''%s''',param.method));
    end
    param.kernel = fftn(ifftshift(param.kernel))/prod(res);
    %fprintf('PSF precomputed with method %s in %f seconds\n',param.method,etime(clock(),t0));
end

data = ifftn(fftn(x,size(param.kernel)).*param.kernel);
y = data(1:res(1),1:res(2));
end
