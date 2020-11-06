function [X, Y] = GenerateScalableSignals(W, Int, Nber, SamplingStrategy)

% Fabien Boux - 11/2020

narginchk(4, 4);

P = length(W);

switch SamplingStrategy
    case 'Grid'
        Y = [];
        
    case 'Random'
        Y  	= Int(1) + (Int(2)-Int(1)) * rand(Nber, P);
        
    case 'qRandom'
        Y 	= Int(1) + (Int(2)-Int(1)) * net(scramble(sobolset(P),'MatousekAffineOwen'), Nber);
        
    otherwise
        error('Invalid sampling strategy')
        
end

X   = ScalableSignals(Y, W);

end