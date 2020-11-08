%% PerformTests.m
%
% Available tests:
%	* 'traj' 		Testing code related to k-space trajectories
%	* 'sens' 		Testing code related to sensitivities
%	* 'nuft accuracy'	Checking accuracy of NUFT implementations
%	* 'nuft speed'		Comparing speed of NUFT implementations
%	* 'noise'		Checking estimation of correlation matrices
%	* 'simu'		Checking consistency of simulation methods
%	* 'dwt'			Checking implementation of the DWT class
%	* 'dct'			Checking implementation of the BlockDCT class
%
%	* 'all'			Run all the above tests
%
% Copyright, Matthieu Guerquin-Kern, 2012

function PerformTests(string)

addpath('test/');
if nargin==0
	string = 'all';
end

switch lower(string)
	case 'all'
		TestGeneratingTraj;pause
		TestGeneratingSens;pause
		TestNUFTaccuracy;pause
		TestNUFTspeed;pause
		TestNoise;pause
		TestSimulation;pause
		TestDWT;pause
		TestBlockDCT;
	case 'traj'
		TestGeneratingTraj;
	case 'sens'
		TestGeneratingSens;
	case 'nuft accuracy'
		TestNUFTaccuracy;
	case 'nuft speed'
		TestNUFTspeed;
	case 'noise'
		TestNoise;
	case 'simu'
		TestSimulation;
	case 'dwt'
		TestDWT;
	case 'dct'
		TestBlockDCT;
	otherwise
		error(sprintf('unknown test ''%s''',string));
end
