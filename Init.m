
%% Description
%
% This script initializes the working environment, ie generates 
% dependencies, creates temporary folders, etc.
% Note that this script is run at the beginning of all the experiments
%
% Fabien Boux - 11/2020


disp(['Running Init.m - ' date])

%% Define path

project_root_path = mfilename('fullpath');
project_root_path = project_root_path(1:end-4);


%% Generate path

addpath(project_root_path)
addpath(genpath([project_root_path filesep 'functions']))
addpath(genpath([project_root_path filesep 'tools']))
addpath(genpath([project_root_path filesep 'Experiments']))


%% Create required folders

if ~exist([project_root_path filesep 'functions'], 'dir'), mkdir([project_root_path filesep 'functions']); end
if ~exist([project_root_path filesep 'temp'], 'dir'), mkdir([project_root_path filesep 'temp']); end 
if ~exist([project_root_path filesep 'figures'], 'dir'), mkdir([project_root_path filesep 'figures']); end 


