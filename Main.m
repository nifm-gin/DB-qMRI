
%% Description
%
% The aim of this script is to run all the others and to produce all
% numerical experiments that allows to illustrate the interest of the
% dictionary-based learning (DBL) method.
%
% Fabien Boux - 01/2020


%% Setting

% Figures are saved in the ./figures folder 
close_figure = 1;


%% Effect of parameter space sampling
% Produce 3 figures:
%   - ParameterSpaceSampling-illustration.fig illustrates
%       2-dimensional projections of 3-dimensional samples obtained from 3
%       different sampling strategies (grid, random and quasi-random)
%
%   - ParameterSpaceSampling.fig
%       Distribution of average RMSE for the different sampling strategies
%       applying the DBL method
%
%   - ParameterSpaceSampling-supp.fig
%       More experiments and average RMSE for both the DBM and DBL methods

clear 
ParameterSpaceSampling
if close_figure == 1, close all; end


%% Impact of the dictionary size
% Produce 2 figures: 
%   - DictionarySize.fig
%       Average RMSE are given as a function of the SNR for different
%       dictionary sizes and both the DBM and DBL methods
%
%   - DictionarySize-supp.fig 
%       More experiments

clear 
DictionarySize
if close_figure == 1, close all; end


%% Boundary behaviour
% Produce 2 figures: 
%   - BoundaryBehaviour.fig
%       Average RMSE over the 2-dimensional space obtained using the DBM
%       and the DBL methods on dictionary smaller than this space.
%
%   - BoundaryBehaviour-supp.fig 
%       Same experiment with few signals added in the dictionary

clear 
BoundaryBehaviour
BoundaryBehaviour_supp
if close_figure == 1, close all; end













