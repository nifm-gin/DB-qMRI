function [X, Y] = GenerateScalableSignals(W, Int, Nber, SamplingStrategy)

% Fabien Boux - 11/2020

narginchk(4, 4);

P   = length(W);

switch SamplingStrategy
    case 'Grid'
        step    = ( Int(2) - Int(1) ) / ( Nber^(1/P) );
        val     = Int(1) + step/2 : step : Int(2) - step/2;
        Y       = arrangement(val,P);    
        
        if Nber ~= length(Y), warning('Not enought signals in grid dictionary'); end                  
                        
    case 'Random'
        Y  	= Int(1) + (Int(2)-Int(1)) * rand(Nber, P);
        
    case 'qRandom'
        Y 	= Int(1) + (Int(2)-Int(1)) * net(scramble(sobolset(P),'MatousekAffineOwen'), Nber);
        
    otherwise
        error('Invalid sampling strategy')        
end

X   = ScalableSignals(Y, W);
        
end