function [qMaps, PD, x, r] = admm_recon(E, data, Dic, ADMM_iter, cg_iter, mu1, mu2, lambda, P, verbose)
% Reconstructs quantitative maps from k-space data by alternately
% solving the inverse imaging problem, constraint to be close to the latest
% dictionary fit, and fitting the series of images to the dictionary.
%
% [qMaps, PD, x]    = admm_recon(E, data, Dic)
% [qMaps, PD, x]    = admm_recon(E, data, Dic, ADMM_iter)
% [qMaps, PD, x]    = admm_recon(E, data, Dic, ADMM_iter, cg_iter)
% [qMaps, PD, x]    = admm_recon(E, data, Dic, ADMM_iter, cg_iter, mu1)
% [qMaps, PD, x]    = admm_recon(E, data, Dic, ADMM_iter, cg_iter, mu1, mu2, lambda)
% [qMaps, PD, x]    = admm_recon(E, data, Dic, ADMM_iter, cg_iter, mu1, mu2, lambda, P)
% [qMaps, PD, x]    = admm_recon(E, data, Dic, ADMM_iter, cg_iter, mu1, mu2, lambda, P, verbose)
% [qMaps, PD, x, r] = admm_recon(___)
%
% Input:
%   E         =  Imaging operator (use LR_nuFFT_operator provided by this
%                toolbox. It can be used for a low rank approximation of the
%                time series, but also for a time frame by time frame
%                reconstruction.
%   data      in [n_samples*nt ncoils]
%                k-space data to be reconstructed. The first dimension
%                represents the readout of all time frames concatted and
%                the second dimension is allows multi-coil data, if E
%                includes the corresponding sensitivity maps.
%   Dic       =  Dictionary struct (see MRF_dictionary.m for details)
%   ADMM_iter =  number of ADMM iterations (default = 10)
%   cg_iter   =  number of conjugate gradient iterations in each ADMM
%                iteration (default = 20)
%   mu1       =  ADMM coupling parameter (dictionary) (default = 1.26e-6,
%                but needs to be changed very likely)
%   mu2       =  ADMM coupling parameter to the spatial regularization
%                term. Has only an effect, if lambda>0 (default = .25)
%   lambda    =  Regularization parameter (default = 0, which results in no
%                spatial regularization)
%   P         =  operator that transforms the images into the space,
%                in which an l21-norm penalty is applied. Has only an
%                effect if lambda>0. Default = 1 (penalty in the image
%                space).
%                Examples:
%                P = wavelet_operator([nx ny nz], 3, 'db2');
%                P = finite_difference_operator([1 2 3]);
%   verbose   =  0 for no output, 1 for plotting the images and in each
%                iteration and give only one output per ADMM iteration in
%                the commandline and 2 for also print the residal of the
%                CG in each CG iteration.
%
%
%
% Output:
%   qMaps = Maps of quantities contained in D.lookup_table
%   PD    = Proton density retrived from the correlation
%   x     = Low rank - or time-series of images
%   r     = residual after all ADMM steps. Use only when you really want to
%           know it since it requires and additional nuFFT operation
%
% For more details, please refer to
%   J. Asslaender, M.A. Cloos, F. Knoll, D.K. Sodickson, J.Hennig and
%   R. Lattanzi, Low Rank Alternating Direction Method of Multipliers
%   Reconstruction for MR Fingerprinting  Magn. Reson. Med., epub
%   ahead of print, 2016.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) Jakob Asslaender, August 2016
% New York University School of Medicine, Center for Biomedical Imaging
% University Medical Center Freiburg, Medical Physics
% jakob.asslaender@nyumc.org
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Manage the input...
% If E is a LR_nuFFT_operator, size(E,2) returns [nx ny (nz) R(nt)]
recon_dim = size(E,2);

if nargin < 4 || isempty(ADMM_iter)
    ADMM_iter = 10;
end
if nargin < 5 || isempty(cg_iter)
    cg_iter = 20;
end
if nargin < 6 || isempty(mu1)
    mu1 = 1.26e-3;
end
if nargin < 7 || isempty(mu2)
    mu2 = .25;
end
if nargin < 8 || isempty(lambda)
    lambda = 0;
end
if nargin < 9 || isempty(P)
    P = 1;
end
if nargin < 10 || isempty(verbose)
    verbose = 1;
end

% ploting stuff
persistent h1 h2 h3 h4 h5
for param = 1:size(Dic.lookup_table,2)
    eval(['persistent h', num2str(param+5)]);
end

%% Initaialize
backprojection = E' * data;
x = 0;
y = zeros(recon_dim);
D = zeros(recon_dim);
if lambda > 0
    Px = zeros(size(P * backprojection));
    z = Px;
    G = Px;
end

r = zeros(1,ADMM_iter);

%% ADMM Iterations
for j=1:ADMM_iter
    tic;
    
    %% update x
    b = backprojection - mu1 * y + mu1 * D .* repmat(sum(conj(D) .* y, length(recon_dim)), [ones(1,length(recon_dim)-1) recon_dim(end)]);
    if lambda > 0
        b = b + mu2 * (P' * (G - z));
        % if j==0, mu1 = 0 in order to realize (DDh)_0 = 1
        f = @(x) E'*(E*x) + (mu1 * (j>1)) * (x - D .* repmat(sum(conj(D) .* x, length(recon_dim)), [ones(1,length(recon_dim)-1) recon_dim(end)])) + (mu2 * (j>1)) * (P' * (P * x));
        %         f = @(x) E'*(E*x) + (mu1 * (j>1)) * (x - D .* repmat(sum(conj(D) .* x, length(recon_dim)), [ones(1,length(recon_dim)-1) recon_dim(end)])) + mu2 * (P' * (P * x));
    else
        % if j==0, mu1 = 0 in order to realize (DDh)_0 = 1
        f = @(x) E'*(E*x) + (mu1 * (j>1)) * (x - D .* repmat(sum(conj(D) .* x, length(recon_dim)), [ones(1,length(recon_dim)-1) recon_dim(end)]));
    end
    x = conjugate_gradient(f,b,1e-6,cg_iter,x,verbose);
    
    
    %% update D and y (and claculate PD and qMaps for the output)
    x  = reshape(x, [prod(recon_dim(1:end-1)), recon_dim(end)]);
    y  = reshape(y, [prod(recon_dim(1:end-1)), recon_dim(end)]);
    D  = reshape(D, [prod(recon_dim(1:end-1)), recon_dim(end)]);
    
    for q=size(x,1):-1:1
        Dx = x(q,:) * Dic.magnetization;
        Dy = y(q,:) * Dic.magnetization;
        [~,idx(q,1)] = max(2*real(Dx.*conj(Dy)) + abs(Dx).^2, [], 2);
    end
    
    DDhx_old = reshape(D .* repmat(sum(conj(D) .* x ,2), [1 recon_dim(end)]), recon_dim);
    
    D = double(Dic.magnetization(:,idx)).';
    Dhx = sum(conj(D) .* x ,2);
    PD = Dhx ./ Dic.normalization(idx).';
    qMaps = Dic.lookup_table(idx,:);
    
    x   = reshape(x,    recon_dim);
    y   = reshape(y,    recon_dim);
    D   = reshape(D,    recon_dim);
    Dhx = reshape(Dhx,  recon_dim(1:end-1));
    PD  = reshape(PD,   recon_dim(1:end-1));
    qMaps = reshape(qMaps, [recon_dim(1:end-1), size(Dic.lookup_table,2)]);
    
    DDhx = D .* repmat(Dhx, [ones(1,length(recon_dim)-1) recon_dim(end)]);
    y = y + x - DDhx;
    
    % Dynamic update of mu1 according to Boyd et al. 2011
    sd = l2_norm(mu1 * (DDhx - DDhx_old));
    rd = l2_norm(x - DDhx);
    %     if rd > 10 * sd
    %         mu1 = 2*mu1;
    %         y = y/2;
    %     elseif sd > 10 * rd
    %         mu1 = mu1/2;
    %         y = 2*y;
    %     end
    
    %% update G and z
    if lambda > 0
        G_old = G;
        Px = P * x;
        G = Px + z;
        Tl2 = l2_norm(G, length(size(G)));
        G = G ./ repmat(Tl2, [ones(1, length(size(G))-1) recon_dim(end)]);
        G = G .* repmat(max(Tl2 - lambda/mu2, 0), [ones(1, length(size(G))-1) recon_dim(end)]);
        G(isnan(G)) = 0;
        z = z + Px - G;
        
        % Dynamic update of mu2 according to Boyd et al. 2011
        rs = l2_norm(Px - G);
        ss = l2_norm(mu2 * (P' * (G - G_old)));
        if rs > 10 * ss
            mu2 = 2*mu2;
            z = z/2;
        elseif ss > 10 * rs
            mu2 = mu2/2;
            z = 2*z;
        end
    end
    
    %% End of actual algorithm
    
    
    %% Caclulate Residuum (only if it is desired)
    if nargout > 3
        r(1,j) = 0.5 * l2_norm(E*DDhx - data).^2;
        if lambda > 0
            r_spatial = sum(abs(col(l2_norm(P * DDhx, length(size(Px))))));
            r(1,j) = r(1,j) + lambda*r_spatial;
        end
    end
    
    
    %% Below here is just plotting stuff...
    if verbose == 1
        % display D*Dh*x and (x-D*Dh*x)
        if (isempty(h1) || ~ishandle(h1)), h1 = figure; end; set(0,'CurrentFigure',h1);
        imagesc34d(abs(    DDhx),0); title([      'D*Dh*x - iteration = ', num2str(j)]); colorbar; colormap gray; axis off;
        if (isempty(h2) || ~ishandle(h2)), h2 = figure; end; set(0,'CurrentFigure',h2);
        imagesc34d(abs(x - DDhx),0); title(['(x - D*Dh*x) - iteration = ', num2str(j)]); colorbar; colormap gray; axis off;
        
        % display PD and T1
        if (isempty(h3) || ~ishandle(h3)), h3 = figure; end; set(0,'CurrentFigure',h3);
        imagesc34d(abs(PD)); colorbar; axis off; title('PD (a.u.)');
        if isfield(Dic, 'plot_details') && length(Dic.plot_details)>size(Dic.lookup_table,2) && ~isempty(Dic.plot_details{size(Dic.lookup_table,2)+1})
            eval(Dic.plot_details{end});
        end
        
        if length(size(PD)) == 2
            for param = 1:size(qMaps,3)
                eval(['if (isempty(h', num2str(param+5), ') || ~ishandle(h', num2str(param+5), ')), h', num2str(param+5), ' = figure; end; set(0,''CurrentFigure'',h', num2str(param+5), ');']);
                imagesc34d(qMaps(:,:  ,param)); colorbar; axis off;
                if isfield(Dic, 'plot_details') && length(Dic.plot_details)>=param && ~isempty(Dic.plot_details{param})
                    eval(Dic.plot_details{param});
                end
            end
        else
            for param = 1:size(qMaps,4)
                eval(['if (isempty(h', num2str(param+5), ') || ~ishandle(h', num2str(param+5), ')), h', num2str(param+5), ' = figure; end; set(0,''CurrentFigure'',h', num2str(param+5), ');']);
                imagesc34d(qMaps(:,:,:,param)); colorbar; axis off;
                if isfield(Dic, 'plot_details') && length(Dic.plot_details)>=param && ~isempty(Dic.plot_details{param})
                    eval(Dic.plot_details{param});
                end
            end
        end
        
        if lambda > 0
            if (isempty(h4) || ~ishandle(h4)), h4 = figure; end; set(0,'CurrentFigure',h4);
            imagesc34d(abs(P'*G), 0, []); title('P''*G'); axis off; colorbar
        end
        
        if nargout > 3
            if (isempty(h5) || ~ishandle(h5)), h5 = figure; end; set(0,'CurrentFigure',h5);
            plot(1:j, log(r(1,1:j)), 'o');
            xlabel('iteration'); ylabel('log(r)'); title('objective function');
        end
        drawnow;
    end
    
    if verbose > 0
        display(['Iteration ', num2str(j)]);
        disp(['primal residual (dictionary) = ', num2str(rd)]);
        disp(['dual   residual (dictionary) = ', num2str(sd)]);
        %         disp(['next mu1 = ', num2str(mu1)]);
        if lambda > 0
            disp(['primal residual (spatial) = ', num2str(rs)]);
            disp(['dual   residual (spatial) = ', num2str(ss)]);
            disp(['next mu2 = ', num2str(mu2)]);
        end
        toc
        disp(' ');
    end
end
end