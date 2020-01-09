
%% Description
%
% Explaination
%
% Fabien Boux - 01/2020


%% Setting

% Execution settings
verbose = 1;
backup  = 1;

% Signal settings

% Experiment settings

% Regression settings
Parameters.K = 50;
Parameters.cstr.Sigma  = 'd*';
Parameters.cstr.Gammat = ''; 
Parameters.cstr.Gammaw = '';
Parameters.Lw = 0;
snr_train  	= inf;


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

