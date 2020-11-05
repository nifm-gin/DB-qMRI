function [bias, var] = BiasVariance(Ytrue, Yestim)


% bias    = nanmean( abs(Ytrue - Yestim) );
% var     = nanmean( ( abs(Ytrue  - Yestim) - nanmean(abs(Ytrue  - Yestim)) ).^2 );


bias    = nanmean( Ytrue - Yestim );
var     = nanmean( ( Yestim - nanmean(Yestim) ).^2 );

end