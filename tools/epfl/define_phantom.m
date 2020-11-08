function [ref, support, mask] = define_phantom(FOV, res)

%simu = 'analytical';
simu = 'rasterized';

%% Defining MR setup

mxsize = res*[1,1];
%FOV = 0.28; % FOV width

%% Define phantom
disp('Defining analytical phantom');
%DefineSimpleBrain;
% DefineBrain; % Analytical computations are pretty long with this one...
%DefineSL;Brain = SL;clear SL;
DefineBrain_fb;
Brain.FOV = FOV*[1, 1];
ref = RasterizePhantom(Brain,mxsize);
mask = (ref>5e-3*max(ref(:)));
mask = imdilate(mask,strel('disk',5));

support = ones(size(ref));