
%% Description
%
% Explaination
%
% Fabien Boux - 01/2020

Init
disp(['Running experiment ' mfilename '.m'])


%% Setting

% Execution settings
verbose = 1;
backup  = 1;

% Signal settings

% Experiment settings

% Regression settings
Model.K = 50;
Model.cstr.Sigma  = 'd*';
Model.cstr.Gammat = ''; 
Model.cstr.Gammaw = '';
Model.Lw = 0;


%% Creating data

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


%% Displaying


%% Exporting figures

if backup == 1
    savefig(fig, ['figures/' mfilename])
end

