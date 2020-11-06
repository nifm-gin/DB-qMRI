function [Estimation, Model] = AnalyzeMRImages(Sequences,Dico,Method,Model,References,outliers,SNRmap)

% Fabien Boux - 11/2020

narginchk(2, 7);
if ~exist('Method','var'),      Method      = 'DB-SL'; end
if ~exist('Model','var'),       Model       = []; end
if ~exist('References','var'),  References  = []; end
if ~exist('SNRmap','var'),      SNRmap      = []; end


% This 'switch' is only required to remain compatible with old notations
% form previous implementations
switch Method
    case 'ClassicMRF'
        Method = 'DBM';
    case 'RegressionMRF'
        Method = 'DB-SL';
    case 'DBL'
        Method = 'DB-SL';
end

% Some formatting
if isempty(Model),	Model = struct(); end 
if iscell(Dico), Dico = Dico{1}; end
if isempty(Sequences), Sequences = Dico.MRSignals; end

% Compute quantification
switch Method
    
    case 'DBM'
        Estimation = qDBM(Sequences, Dico, References);
        
    case 'DB-SL'
        [Estimation, Model] = qDBSL(Sequences, Dico, References, Model, SNRmap);
        
    case 'DB-DL'
        [Estimation, Model] = qDBDL(Sequences, Dico, References, Model);
end

end
        