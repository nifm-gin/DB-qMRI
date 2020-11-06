function [Theta] = UpdateSigma(Theta, var_noise)

narginchk(2, 2);

Theta.Sigma = Theta.Sigma + var_noise * eye( size(Theta.Sigma(:,:,1)) );

end