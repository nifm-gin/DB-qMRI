%% Sinusoidal2DMatrix.m
%
% Function that provides the matrix that links the vector of parameters s
% to the values of a weighted sum of low-frequency complex exponentials.
% The corresponding frequencies, localized on a grid that is sampled twice
% finer than the Shannon grid (sampling step is the inverse of the FOV),
% define a centered square whose width is parametrized by L.
%
% INPUT:    * R: (2,m) matrix that defines m 2D points of evaluation
%               (coordinates in grid units, i.e. between -0.5 and 0.5)
%           * L: width of the square Fourier support
%
% OUTPUT:   * M: (m,L*L) matrix that links the vector of parameters s 
%               defining the linear combination of sinusoidals, to the
%               values it takes at the m points.
%
% Matthieu Guerquin-Kern, Biomedical Imaging Group / EPF Lausanne,
%  spring 2011
function M = Sinusoidal2DMatrix(R,L,SPEEDUP)
if nargin<3
    SPEEDUP = true;
end

[kx,ky] = GenerateFullCart2DGrid(L*[1,1]);
if SPEEDUP
    ex = exp(1j*pi*R(1,:)).';
    ey = exp(1j*pi*R(2,:)).';
    vx = cell(1,floor(L/2)+1);
    vy = cell(1,floor(L/2)+1);
    vx{1} = ones(size(R,2),1);
    vy{1} = ones(size(R,2),1);
    for i = 2:floor(L/2)+1
        vx{i} = vx{i-1}.*ex;
        vy{i} = vy{i-1}.*ey;
    end
    clear ex ey;
    M = ones([size(R,2),numel(kx)]);
    for i = 1:numel(kx)
        if kx(i)<0
            M(:,i) = conj(vx{abs(kx(i))+1});
        else
            M(:,i) = vx{abs(kx(i))+1};
        end
        if ky(i)<0
            M(:,i) = M(:,i).*conj(vy{abs(ky(i))+1});
        else
            M(:,i) = M(:,i).*vy{abs(ky(i))+1};
        end
    end
else
    M = exp(1j*pi*(kx(:)*R(1,:)+ky(:)*R(2,:))).';
end
