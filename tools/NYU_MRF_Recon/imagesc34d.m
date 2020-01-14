function im = imagesc34d(X, rescale, rc, varargin)
% Function to 2D, 3D or 4D data
%
% imagesc34d(X)
% imagesc34d(X, rescale)
% imagesc34d(X, rescale, rc)
% imagesc34d(X, rescale, rc, varargin)
%
%
% Input:
%   X        = Array of data that can have 2-4 dimension:
%              2D data: This function is equivalent to imagesc
%              3D data: A mosaic is formed with approx. as many rows as
%                       columns, if rc not specified.
%              4D data: A mosaic is formed where dim. 3 is plotted in
%                       columns and dim. 4 is plotted in rows, if rc not
%                       specified.
%   rescale  = 0: Scale all mosaic images the same
%              1: Devide all mosaic images by its maximum
%   rc       = [rows colums]: Specify the number of rows and columns
%   varargin = All inputs exceeding 3 are forwarded to imagesc
%
% If X is complex valued, its absolute value is plotted instead.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) Jakob Asslaender, August 2016
% New York University School of Medicine, Center for Biomedical Imaging
% University Medical Center Freiburg, Medical Physics
% jakob.asslaender@nyumc.org
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


nr = size(X,1);
nc = size(X,2);
slices = size(X,3)*size(X,4);

if nargin<2 || isempty(rescale)
    rescale=0;
end

if nargin<3 || isempty(rc)
    if length(size(X))==4
        cols = size(X,3);
        rows = size(X,4);
    else
        rows = floor(sqrt(slices));
        cols = ceil(slices / rows);
    end
elseif prod(rc)<slices
    error('Number of rows and columns is too small for specified data.');
else
    rows = rc(1);
    cols = rc(2);
end


if length(size(X))==4
    X = X(:,:,:);
end

M = zeros(rows*nr, cols*nc);

c = 1;
r = 1;
for k=1:slices
    if rescale
        M((r-1)*nr+1:r*nr, (c-1)*nc+1:c*nc) = X(:,:,k)/max(col(X(:,:,k)));
    else
        M((r-1)*nr+1:r*nr, (c-1)*nc+1:c*nc) = X(:,:,k);
    end
    c = c + 1;
    if c > cols
        c = 1;
        r = r + 1;
    end
end

if ~isreal(M)
    M = abs(M);
end

if nargout == 0
    if isempty(varargin)
        imagesc(M);
    else
        imagesc(M, varargin{:});
    end
else
    if isempty(varargin)
        im = imagesc(M);
    else
        im = imagesc(M, varargin{:});
    end
end

end
