function x = ifft_nd(x,dim)
% Performes an inverse fast Fourier transformation of x along the dimensions 
% specified in dim. In agreement with the usual definition of k-space in
% MRI, fftshifts are applied before and after the FFT.
%
% x = ifft_nd(x)
% x = ifft_nd(x,dim)
%
% Input:
%   x   = multi dimensional array of data
%   dim = vector of dimensions along which the FFT is performed
%                  (optional; default = all dimensions)
%
% Example:
%   img  = repmat(phantom(256), [1 1 4]); % 4 repetitions of shepp logan 
%   k    =  fft_nd(img, [1 2]);            % k-space of each of them
%   img2 = ifft_nd(img, [1 2]);            % the original image
%
% see also fft_nd
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (c) Jakob Asslaender, August 2016
% New York University School of Medicine, Center for Biomedical Imaging
% University Medical Center Freiburg, Medical Physics
% jakob.asslaender@nyumc.org
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 2 || isempty(dim)
    dim = 1:length(size(x));
end

for j=1:length(dim)
    x=fftshift(ifft(fftshift(x,dim(j)),[],dim(j)),dim(j));
end