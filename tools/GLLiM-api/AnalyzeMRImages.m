function [Estimation, Model] = AnalyzeMRImages(Sequences,Dico,Method,Model,References,outliers,SNRmap)

if nargin < 3, error('Not enought input arguments'); end
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
if isempty(Sequences), Sequences = Dico{1}.MRSignals; end
if isempty(Model),	Model = struct(); end 
if iscell(Dico), Dico = Dico{1}; end


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
        