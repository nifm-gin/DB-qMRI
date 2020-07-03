function [] = OverlayImages(anatImg, paramMap, roiMask, colorMap, bounds)

if ~exist('colorMap','var'),    colorMap = colormap; end
if ~exist('bounds','var'),      bounds   = []; end


paramMap(~roiMask) = nan;

anatImg     = anatImg/(max(anatImg(:)));  	% normalize base (anatomical) image
rgbSlice    = anatImg(:,:,[1 1 1]);      	% converting to RGB (ignore colormaps)

imshow(paramMap, 'Colormap', colorMap);   	% show parametric image
colorbar;
if ~isempty(bounds)
    caxis(bounds)
else
    caxis([min(paramMap(:)) max(paramMap(:))])
end
hold on;

h           = imshow(rgbSlice); % superimpose anatomical image
set(h, 'AlphaData', double(~roiMask));              % make pixels in the ROI transparent
  