function [Yestim, Cov, Kurt, Pik, Mahaldist] = EstimateParametersFromModel(X, f_estim, verb)

narginchk(2, 3);
if nargin == 2, verb = 0; end

[Yestim, Cov, Pik, Mahaldist] = gllim_inverse_map_and_cov(X', f_estim, verb);

Kurt = [];
% [Yestim, Cov, Kurt] = gllim_inverse_map_and_kurt(X', f_estim, verb);

Yestim = Yestim';