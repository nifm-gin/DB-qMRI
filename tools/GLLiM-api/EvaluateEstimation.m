function [Rmse, Nrmse, Mae, Nmae] = EvaluateEstimation(Ytrue, Yestim)

Rmse    = nanmean( (Ytrue - Yestim).^2 ).^.5;
Nrmse   = Rmse ./ Rmse;

Mae     = nanmean( abs(Ytrue - Yestim) );
Nmae    = Mae ./ Mae;

end