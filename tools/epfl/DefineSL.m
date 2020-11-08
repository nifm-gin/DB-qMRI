%% DefineSL.m
%
% SEE: phantom
%
% Matthieu Guerquin-Kern, Biomedical Imaging Group / EPF Lausanne,
% 30-10-2009 (dd-mm-yyyy)

SL.FOV = 1*[1,1]; % the FOV actually changes the shape of the phantom
%   This head phantom is the same as the Shepp-Logan except
%   the intensities are changed to yield higher contrast in
%   the image.  Taken from Toft, 199-200.
%
%         A    a     b    x0    y0    phi
%        ---------------------------------
shep = [  1   .69   .92    0     0     0
        -.8  .6624 .8740   0  -.0184   0
        -.2  .1100 .3100  .22    0    -18
        -.2  .1600 .4100 -.22    0     18
         .1  .2100 .2500   0    .35    0
         .1  .0460 .0460   0    .1     0
         .1  .0460 .0460   0   -.1     0
         .1  .0460 .0230 -.08  -.605   0
         .1  .0230 .0230   0   -.606   0
         .1  .0230 .0460  .06  -.605   0   ];

SL.region = cell(1,size(shep,1));
for i = 1:size(shep,1)
    SL.region{i} = struct('type','ellipse',...
        'center',[-shep(i,5),shep(i,4)]/2,... % in FOV units
        'angle',shep(i,end)*pi/180,...
        'weight',shep(i,1),...
        'width',[shep(i,[3,2])]); % in FOV units
end
clear shep i;