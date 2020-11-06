function [Rmse, Nrmse, Mae, Nmae] = EvaluateEstimation(Ytrue, Yestim)

% Compute accuracy metrics
%
% Fabien Boux - 10/2018

narginchk(2, 2);

Rmse    = nanmean( (Ytrue - Yestim).^2 ).^.5;
Nrmse   = Rmse ./  nanmean( ( Ytrue - nanmean(Ytrue)).^2).^.5;

Mae     = nanmean( abs(Ytrue - Yestim) );
Nmae    = Mae ./ nanmean( Ytrue - nanmean(Ytrue) );

end