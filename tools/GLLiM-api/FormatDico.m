function [Dico] = FormatDico(Xtrain, Ytrain, Tacq)

% Fabien Boux - 11/2020

narginchk(2, 3);

if size(Xtrain,1) ~= size(Ytrain, 1)
    error('Sizes do not match')
end

Dico.MRSignals = Xtrain;
Dico.Parameters.Par = Ytrain;

if exist('Tacq','var')
    Dico.Tacq = Tacq;
end