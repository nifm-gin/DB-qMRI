function [h_hat, g_hat, htld_hat, gtld_hat] = BiorthWavFiltersSep(N, J, decimation, family)
% [h_hat, g_hat, htld_hat, gtld_hat] = BiorthWavFiltersSep(N, J, decimation, family)
%
% Generate the DFTs of the scaling and wavelet filters for various separable
% wavelet families, as required by the decomposition and reconstruction routines.
%
% The following features are supported:
% - arbitrary dimensionality;
% - orthonormal or biorthogonal bases;
% - dyadically decimated or undecimated representations;
% - possibly non-stationary refinement filters;
% - dimension-dependent wavelet families and decomposition depths.
%
% N:          Dimension-specific filter length.
% J:          Dimension-specific decomposition depth.
% decimation: True or false.
% family:     Wavelet family (see BiorthWavFilters1D() routine or,
%             for specifying dimension-dependent families, see code).
%
% h_hat:      DFT of the synthesis scaling filter.
% g_hat:      DFT of the synthesis wavelet filter.
% htld_hat:   DFT of the analysis scaling filter.
% gtld_hat:   DFT of the analysis wavelet filter.
%
% NOTE: the output variables are 2D cell arrays whose dimensions are:
%       - the largest decomposition depth;
%       - the dimensionality of N.
%
% (c) Cedric Vonesch, 2006.10.06-2008.04.04

% Decomposition depth
Jmax = max(J);

% Dimensionality
D = length(N);

% Dimension-dependent wavelet families
wavelet = cell(D, 1);
switch lower(family)
    case 'spline4haar'
        for d = 1:D-1
            wavelet{d} = 'Spline4';
        end
        wavelet{D} = 'Haar';
    case 'spline2haar'
        for d = 1:D-1
            wavelet{d} = 'Spline2';
        end
        wavelet{D} = 'Haar';
    case 'spline4spline3'
        for d = 1:D-1
            wavelet{d} = 'Spline4';
        end
        wavelet{D} = 'Spline3';
    otherwise
        for d = 1:D
            wavelet{d} = family;
        end
end

% Filters
h_hat = cell(Jmax, D);
g_hat = cell(Jmax, D);
htld_hat = cell(Jmax, D);
gtld_hat = cell(Jmax, D);
for j = 1:Jmax
    for d = 1:D
        if j <= J(d)
            [h_hat{j, d}, g_hat{j, d}, htld_hat{j, d}, gtld_hat{j, d}] = BiorthWavFilters1D(N(d)/2^(j-1), wavelet{d}, j);
            if ~decimation(1)
                h_hat{j, d} = h_hat{j, d} / 2; % For synthesis (direct)
                h_hat{j, d} = repmat(h_hat{j, d}, 2^(j-1), 1);
                htld_hat{j, d} = repmat(htld_hat{j, d}, 2^(j-1), 1); % For synthesis (adjoint)
                if decimation(2)
                    disp('youyou')
                    g_hat{j, d} = g_hat{j, d} / 2; % For synthesis (direct)
                    
                end
                g_hat{j, d} = repmat(g_hat{j, d}, 2^(j-1), 1);
                gtld_hat{j, d} = repmat(gtld_hat{j, d}, 2^(j-1), 1); % For synthesis (adjoint)
            end
        end
    end
end