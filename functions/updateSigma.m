function [Theta_updated] = updateSigma(Theta, var_noise)

Theta_updated = Theta;
Theta_updated.Sigma = Theta_updated.Sigma + var_noise * eye(size(Theta_updated.Sigma(:,:,1)));
end

