%% TestGeneratingTraj.m
%
% Generating the different kinds of trajectories and scaling them to grid
% units for visualization.
%
% Copyright, Matthieu Guerquin-Kern, 2012

%clear all;clc;
close all;
disp('Performing the test for k-space trajectories...');
FOV = .24;
mxsize = 256;
mx = mxsize*[1,1];

%% Cartesian
w = GenerateCartesianTraj(FOV, FOV/mxsize, 1, 1);
[k,mx] = TrajInGridUnits(w,FOV);
figure(99);plot(k(:,1),k(:,2),'.-');axis square;title(sprintf('Cartesian traj: %d x %d',mx(1),mx(2)));

%% Spiral
w=GenerateSpiralTraj(FOV,FOV/mxsize,1,1,10,1, .031, 200,true,false);
[k,mx] = TrajInGridUnits(w,FOV);
figure(98);plot(k(:,1),k(:,2),'.-');axis square;title(sprintf('Spiral traj: %d x %d',mx(1),mx(2)));

%% Radial
w=GenerateRadialTraj(FOV,FOV/mxsize,1,1);
[k,mx] = TrajInGridUnits(w,FOV);
figure(97);plot(k(:,1),k(:,2),'.');axis square;title(sprintf('Radial traj: %d x %d',mx(1),mx(2)));
