
%% Description
% Explaination
%
% Fabien Boux 01/2020


%% Setting

% Execution settings
verbose = 1;
backup  = 1;

% Signal settings

% Experiment settings

% Regression settings


%% Data Creation

% Adding to path
addpath(genpath('functions'))
addpath(genpath('tools'))

% Creating 


%% Processing 


%% Saving 

if backup == 1
    clear tmp*
    save(['temp/' mfilename])
end


%% Display


%% Exporting figures

if backup == 1
    savefig(fig, ['outputs/' mfilename])
end

