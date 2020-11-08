%% TestSimulation.m
%
% In a typical single channel MRI setting, the consistency of two different
% MRI data simulation methods is checked in both k-space and reconstructed image
% domain.
%
% Copyright, Matthieu Guerquin-Kern, 2012

close all;
disp('Performing the test for the MR data synthesis methods...');
% simulation parameters
DefineSimpleBrain;
mxsize = 64*[1,1];
R = 1;
f_sampling = 1/1.3;
Ncoils = 1;
FOV = 0.28; % FOV width
Brain.FOV = FOV*[1, 1];
pixelsize = FOV./mxsize;
w = GenerateSpiralTraj( FOV, pixelsize, f_sampling, R, 10, 1, 0.031, 200, false, false);
%w = GenerateCartesianTraj(FOV, pixelsize, f_sampling, R);

sensitivity = GenerateSensitivityMap( FOV, pixelsize, Ncoils, .07, .17);
im = RasterizePhantom(Brain,mxsize);
kernel = ones(7);kernel = kernel/sum(kernel(:));
support = (conv2(im,kernel,'same')>5e-3*max(im(:)));

%%
disp('Simulating analytical MR data');
m_ana = zeros(size(w,1),Ncoils);
sens=cell(1,Ncoils);
for indCoil = 1:Ncoils
    sens{indCoil} = SensFitting(sensitivity(:,:,indCoil),'sinusoidal',5,support);
    m_ana(:,indCoil) = MRDataAnalytical(Brain, sens{indCoil}, w);
end
m_ana = prod(mxsize).*m_ana;
%%
disp('Simulating rasterized MR data');
refinement = 8; % to simulate a continuous setting with rasterized data
sens_rast = zeros([refinement*mxsize Ncoils]);
im_rast = RasterizePhantom(Brain,refinement*mxsize);
m_rast = zeros(size(w,1),Ncoils);
for indCoil = 1:Ncoils
    [im_raster,sens_rast(:,:,indCoil)] = RasterizePhantom(Brain,refinement*mxsize,sens{indCoil});
	m_rast(:,indCoil) = MRDataRasterized(sens_rast(:,:,indCoil).*im_rast, w, FOV)/refinement^ndims(im_rast);
end

%% Comparison
fprintf('the following number should be as close to 0 as possible: %g\n',norm(m_rast(:)-m_ana(:))/norm(m_ana(:)));

%% Recons
k = TrajInGridUnits(w, FOV, mxsize);

[a,A,P] = Prepare4Recon(m_ana, k, sensitivity, support);
x_ana = ReconCG(a,@(x) A(x) + 0*x,zeros(mxsize),60,P);
figure(5);imagesc((abs(x_ana)).^.8);colormap gray;axis image;colorbar;title('reconstructed image (analytical)');

[a,A,P] = Prepare4Recon(m_rast, k, sensitivity, support);
x_rast = ReconCG(a,@(x) A(x) + 0*x,zeros(mxsize),60,P);
figure(6);imagesc((abs(x_rast)).^.8);colormap gray;axis image;colorbar;title('reconstructed image (rasterized)');
figure(7);imagesc((abs(x_ana-x_rast)).^.8);colormap(1-gray);axis image;colorbar;title('reconstruction error with inverted gray levels');