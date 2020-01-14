classdef wavelet_operator
    % wavelet_operator immitates a matrix that performes a wavelet
    % transformation.
    %   The constuctor takes the size of the dimensions along which the
    %   wavelet transformation shall be perforemed (e.g. [nx ny nz]). The
    %   wavelet transformation is performed when multiplying the operator
    %   (mtimes) to an image. The image can also have one dimension more
    %   than the operator (e.g. time). In this case, the wavelet
    %   transformation is performed for each frame.
    %   E.g. in order to calculate a gradient, the adjoint of the operator
    %   can be use by calling A'.
    %
    %   Example: A = wavelet_operator([nx ny], 3, 'db2');
    %            x = rand(256, 256, 10); % the last dimension is time
    %            d = A*x;                % apply WL along dim 1 and 2
    %            g = A'*d;               % adjoint operation
    %
    % see also wavedec, wavedec2, wavedec3
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % (c) Jakob Asslaender, August 2016
    % New York University School of Medicine, Center for Biomedical Imaging
    % University Medical Center Freiburg, Medical Physics
    % jakob.asslaender@nyumc.org
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    properties
        space = 0;      % Dimensionality of the wavelet transformation
        params = [];    % Bunch of wavelet parameters
        dim = [];       % Size of each image to be transformed
        sizes_all = []; % For the 3D wavelet
        N = [];         % Order of wavelet transformation
        wname = '';     % Name of the wavelet
        adjoint = 0;    % Boolean indicating if the matrix is adjoint
        rand_shift = [];
    end
    
    methods
        function W = wavelet_operator(dim, N, wname)
            % Constructor:
            % W = waveletDecompositionOperator(dim, N, wname)
            % Input:   dim = dimension of image (e.g. [64 64 32])
            %            N = order of decomposition
            %        wname = name of wavelet
            % Output: Object
            %
            % see also wavedec, wavedec2, wavedec3
            
            if nargin==0
                W.space = 0;
            else
                if length(dim)==1
                    W.space = 1;
                elseif length(dim)==2
                    if dim(1)==1 || dim(2)==1
                        W.space = 1;
                    else
                        W.space = 2;
                    end
                elseif length(dim)==3
                    W.space = 3;
                end
            end
            
            W.rand_shift = zeros(1,length(dim));
            
            if W.space==0
                W.params.sizes = [];
                W.params.sizeINI = [];
                W.params.level = [];
                W.params.mode = [];
                W.params.filters = [];
                W.dim = [];
                W.space = [];
                W.sizes_all = [];
            elseif W.space==1
                [~,sizes] = wavedec(zeros(dim,1),N,wname);
                W.params.sizes = sizes;
                W.params.sizeINI = [];
                W.params.level = [];
                W.params.mode = [];
                W.params.filters = [];
                W.dim = [dim 1];
                W.space = 1;
                W.sizes_all = [];
            elseif W.space==2
                [~,sizes] = wavedec2(zeros(dim),N,wname);
                W.params.sizes = sizes;
                W.params.sizeINI = [];
                W.params.level = [];
                W.params.mode = [];
                W.params.filters = [];
                W.dim = dim;
                W.space = 2;
                W.sizes_all = [];
            elseif W.space==3
                W.params = wavedec3(zeros(dim),N,wname);
                W.dim = dim;
                W.space = 3;
                sizes_all = zeros(length(W.params.dec),3);
                for k=1:length(W.params.dec)
                    sizes_all(k,:) = size(W.params.dec{k});
                end
                W.sizes_all = sizes_all;
            else
                error('dimension is not specified properly.');
            end
            if nargin==0
                W.N = [];
                W.wname = '';
                W.adjoint = 0;
            else
                W.N = N;
                W.wname = wname;
                W.adjoint = 0;
            end
        end
        
        function W = ctranspose(W)
            % Sets the objection in the adjoint mode.
            % Called when using A'
            W.adjoint = ~W.adjoint;
        end
        
        function Q = mtimes(W,B)
            % Performs the actual wavelet calculation.
            %   Called when using A*B or A'*B
            %   Input: A is usually the wavelet operator
            %          B is the image (adjoint=0) or a vector of wavelet
            %            components (adjoint=1)
            
            if isa(W, 'wavelet_operator')
                if W.adjoint
                    B = circshift(B, W.rand_shift);
                    if W.space==1
                        for i=size(B,2):-1:1
                            Q(:,i) = conj(waverec(B(:,i),W.params.sizes,W.wname));
                        end
                    elseif W.space==2
                        for i=size(B,2):-1:1
                            Q(:,:,i) = conj(waverec2(B(:,i),W.params.sizes,W.wname));
                        end
                    elseif W.space==3
                        X.sizeINI = W.params.sizeINI;
                        X.level = W.params.level;
                        X.mode = W.params.mode;
                        X.filters = W.params.filters;
                        X.sizes = W.params.sizes;
                        ps = prod(W.sizes_all,2);
                        Y = cell(1,length(ps));
                        for i=size(B,2):-1:1
                            startpt = 1;
                            for k=1:length(ps)
                                Y{k} = reshape(B(startpt:startpt+ps(k)-1,i),W.sizes_all(k,:));
                                startpt = startpt + ps(k);
                            end
                            X.dec = Y;
                            Q(:,:,:,i) = waverec3(X);
                        end
                    end
                else
                    if W.space==1
                        for i=size(B,2):-1:1
                            Q(:,i) = wavedec(B(:,i),W.N,W.wname);
                        end
                    elseif W.space==2
                        for i=size(B,3):-1:1
                            Q(:,i) = wavedec2(B(:,:,i),W.N,W.wname).';
                        end
                    elseif W.space==3
                        for i=size(B,4):-1:1
                            X = wavedec3(B(:,:,:,i),W.params.level,W.wname);
                            X = X.dec;
                            if i == size(B,4)
                                ps = prod(W.sizes_all,2);
                                Q = zeros(sum(ps),size(B,4));
                            end
                            startpt = 1;
                            for k=1:length(X)
                                Q(startpt:startpt+ps(k)-1,i) = col(X{k});
                                startpt = startpt + ps(k);
                            end
                        end
                    end
                    Q = circshift(Q, -W.rand_shift);
                end
            else % now B is the operator and A is the vector
                Q = conj(mtimes(B',conj(W)));
                
            end
        end
        
        function s = size(W,n)
            if nargin < 2
                s(2) = prod(W.dim);
                s(1) = sum(prod(W.sizes_all,2));
            else
                if n==1
                    s = sum(prod(W.sizes_all,2));
                elseif n==2
                    s = W.dim;
                else
                    s = 1;
                end
            end
        end
        
        function A = update_rand_shift(A)
            A.rand_shift = round(rand(1,length(A.dim)).*A.dim);
        end
    end
end