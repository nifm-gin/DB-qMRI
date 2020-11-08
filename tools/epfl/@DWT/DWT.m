% Class of object that represents discrete wavelet transform. The classical
% operators for matrix multiplication are overloaded to permit some
% manipulations as for matrices. The wavelet coefficients are represented
% in the format WAVELET.
%
% The class relies on the wavelet implementations of Cedric Vonesch in the
% Biomedical Imaging Group. In particular, the operations are performed in
% the fourier domain.
% Cf. http://bigwww.epfl.ch/algorithms/mltldeconvolution/
%
% See also: WAVELET
%
% Matthieu Guerquin-Kern, Biomedical Imaging Group / EPF Lausanne,
% 17-07-2009 (dd-mm-yyyy)

classdef DWT
    properties
        adjoint = 0;
        inverse = 0;
        J = 5;
        family = 'haar'; % List: of families: cf BiorthWavFilters1D.m
        decimation = true; % true or false
        size; % size of the object to be transformed. The number of dimensions: numel(size)
        h_hat;
        g_hat;
        htld_hat;
        gtld_hat;
        hstr_hat;
        gstr_hat;
        htldstr_hat;
        gtldstr_hat;
    end
    
    methods
        function  A = DWT(J,res, decimation, family)    % constructor
            A.adjoint = 0;
            A.inverse = 0;
            ind = find(res>1);
            A.size = res(ind);
            if numel(J)<numel(res)
                J = J(1)*ones(1,numel(res));
            end
            A.J = J;
            switch length(decimation)
                case 1
                    A.decimation = [decimation,decimation];
                case 2
                    A.decimation = decimation;
                otherwise
                    error('DWT: the decimation parameter should have length 1 or 2');
            end
            
            
            A.family = family;
            
            switch A.family
                %case 'haar'
                %    [A.h_hat,A.g_hat,A.htld_hat,A.gtld_hat] = wfilters('haar');
                otherwise
                    % Precomputation of the filters
                    [A.h_hat, A.g_hat, A.htld_hat, A.gtld_hat] = BiorthWavFiltersSep(res, J, decimation, family);
                    A.hstr_hat = cellfuncell(@conj, A.h_hat); % For analysis (adjoint)
                    A.gstr_hat = cellfuncell(@conj, A.g_hat); % For analysis (adjoint)
                    A.htldstr_hat = cellfuncell(@conj, A.htld_hat); % For analysis (direct)
                    A.gtldstr_hat = cellfuncell(@conj, A.gtld_hat); % For analysis (direct)
            end
        end
        
        function s = display(A)                             % description display
            if A.inverse ==0 % DWT
                if A.adjoint == 0 %A*x
                    direction = 'direct';
                    image = 'input';
                    filter = [];
                else %At*x
                    direction = 'inverse';
                    image = 'output';
                    filter = ' using dual filters';
                end
            else % DWT inv
                if A.adjoint == 0 %A^-1*x
                    direction = 'inverse';
                    image = 'output';
                    filter = [];
                else %A^-t*x
                    direction = 'direct';
                    image = 'input';
                    filter = ' using dual filters';
                end
            end
            switch A.decimation(1)
                case 1, decim1_strg = '';
                case 0, decim1_strg = ', undecimated coarse coeff,';
            end
            switch A.decimation(2)
                case 1, decim2_strg = '';
                case 0, decim2_strg = ', undecimated wavelet coeff,';
            end
            size_strg = [];
            for i = 1:length(A.size)
                size_strg = [size_strg 'x' num2str(A.size(i))];
            end
            depth_strg = [];
            for i = 1:length(A.J)
                depth_strg = [depth_strg 'x' num2str(A.J(i))];
            end
            s =[direction ' ' num2str(length(A.size)) '-D ' A.family decim1_strg decim2_strg ' DWT ' 'of depth ' depth_strg(2:end) ' for an ' image ' of size ' size_strg(2:end) filter ];
            disp(s);
        end
    end % methods
end % classdef