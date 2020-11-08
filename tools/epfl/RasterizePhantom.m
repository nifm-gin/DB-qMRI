%% RasterizePhantom.m
%
% Function that rasterizes a 2D analytical phantom composed of ellipses,
% polygons and Bezier-defined regions. Supports a sensitivity weighting.
% If no ouput is required, displays the rasterized phantom.
%
% INPUTS:   * phantom object
%           * resolution for rasterization
%           * sensitivity (optional)
%
% OUTPUTS:  * rasterized image
%           * rasterized sensitivity
%
% Matthieu Guerquin-Kern, Biomedical Imaging Group / EPF Lausanne,
% 20-10-2009 (dd-mm-yyyy)
% revised in Feb. 2011

function [im,sens_rast] = RasterizePhantom(phantom,res,s,MATLAB_SL_CONSISTENT)
%Input check
if nargin<4, MATLAB_SL_CONSISTENT = false;end
if nargin<3, s = 1;end
if nargin<2, res = floor(256);end
if numel(res)==1, res = [res,res];end

% Use Bruno Luong's insidepoly.m
% (see
% http://www.mathworks.com/matlabcentral/fileexchange/27840-2d-polygon-inte
% rior-detection)
% if possible or fallback to matlab's inpolygon
if exist('insidepoly_dblengine','file')==3, inpoly = @insidepoly;
else inpoly = @inpolygon;end

% Rasterization of the phantom
im = zeros(res);
[X,Y] = GenerateFullCart2DGrid(res);

if MATLAB_SL_CONSISTENT % this spatial sampling allows consistency with Matlab's internal Shepp-Logan but the origin is not sampled for even resolutions.
    disp('using axis in consistency with the Matlab phantom.m')
    X = (X+.5*mod(res(1)+1,2))/(res(1)-1);
    Y = (Y+.5*mod(res(2)+1,2))/(res(2)-1);
else
    X = X/res(1);
    Y = Y/res(2);
end

Nreg = length(phantom.region);
use_waitbar = false;%Nreg>5;
if use_waitbar
    hbar = waitbar(0,'RASTERIZATION');
end
for indRegion = 1:Nreg
    if use_waitbar
        waitbar((indRegion-1)/Nreg,hbar); %,sprintf('RASTERIZATION -> computing region %d of %d',indRegion,Nreg)
    end
    region = phantom.region{indRegion};
    switch region.type
        case {'ellipse'}
            ct = cos(region.angle);
            st = sin(region.angle);
            x1 = X-region.center(1);
            x2 = Y-region.center(2);
            u1 = 2/(region.width(1))*(ct*x1+st*x2);
            u2 = 2/(region.width(2))*(-st*x1+ct*x2);
            support_ell = sqrt(u1.^2+u2.^2)<=1;
            im(support_ell) = im(support_ell)+region.weight;
        case{'polygon'}
            im = im + region.weight*double(inpoly(X,Y,region.vertex(:,1),region.vertex(:,2)));
        case {'bezier'}
            im = im + region.weight*double(InBezierDefinedRegion(region,X,Y));
        otherwise
            error('''%s'' is an unsupported type of region',region.type);
    end
end
if use_waitbar
    close(hbar);
end

% Sensitivity weighting
if isnumeric(s)&&numel(s)==1
    sens_rast = s*ones(res);
elseif isnumeric(s)&&size(s)==res
    sens_rast = s;
elseif isstruct(s)
    sens_rast = zeros(res);
    [X1,X2] = GenerateFullCart2DGrid(res);
    test = 0;
    R = [(X1(:)'+test*.5)/(res(1)-test*1);(X2(:)'+test*.5)/(res(2)-test*1)];clear X1 X2;
    switch s.model
        case 'polynomial'
            D = s.param;%floor((sqrt(1+numel(s)*8)-3)/2); % works only in 2D!
            M = Polynomial2DMatrix(R,D);
            sens_rast(:) = M*s.data;clear M;
        case 'sinusoidal'
            L = s.param;
            M = Sinusoidal2DMatrix(R,L);clear R;
            sens_rast(:) = M*s.data;clear M;
    end
else
    error('unsupported type of sensitivity');
end
im = im.*sens_rast;

% Output Check
if nargout==0
    imagesc(abs(im));colormap gray;colorbar;axis image;drawnow;
    title(sprintf('%dx%d rasterization of the phantom',res(1),res(2)));
    x = (-0.4:0.2:0.4);
    set(gca,'XTick',(x+0.5)*res(1),'XTickLabel',x*phantom.FOV(1),'YTick',(x+0.5)*res(2),'YTickLabel',x*phantom.FOV(2),'XDir','normal','YDir','reverse');
end

%% Subfunctions
    function map = InBezierDefinedRegion(region,X,Y)
        if exist('insidepoly_dblengine','file')==3, inpoly = @insidepoly;
	else inpoly = @inpolygon;end
        r = control2node(region.control);
        map = double(inpoly(X,Y,r(:,1),r(:,2))); % Fill the polygon of vertices rn
        c = circshift(region.control,[1,0]);
        rp1 = circshift(r,[-1,0]);
        beta = c-circshift(c,[1,0]);
        gamma = rp1+r-2*c;
        a = beta(:,1).*gamma(:,2)-beta(:,2).*gamma(:,1);
        for i = 1:length(a)
            ind = find(inpoly(X,Y,[r(i,1) c(i,1) rp1(i,1)],[r(i,2) c(i,2) rp1(i,2)])); % consider the points inside the triangle
            b = -(X(ind)-r(i,1))*gamma(i,2)+(Y(ind)-r(i,2))*gamma(i,1);
            d = -(X(ind)-r(i,1))*beta(i,2)+(Y(ind)-r(i,2))*beta(i,1);
            map(ind(b.^2<a(i)*d)) = (a(i)>=0);
            % a>=0 for outward-pointing triangles: add the interior points
            % a<0 for inside-pointing triangles: remove the exterior points
        end
    end

end
