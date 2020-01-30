function [Yestim, score, Xestim] = EstimateParametersFromGrid(Xobs, Xgrid, Ygrid, verb)

narginchk(3, 4)
if nargin == 3, verb = 0; end

% normalization
Xobs_normalized     = (1 ./ sum(Xobs.^2, 2) .^.5) * ones(1, size(Xobs, 2))  .* Xobs;
Xgrid_normalized    = (1 ./ sum(Xgrid.^2, 2).^.5) * ones(1, size(Xgrid, 2)) .* Xgrid;

% dot-product/scalar product comparison
score 	= Xobs_normalized * conj(Xgrid_normalized)';
[~, idx] = max(score, [], 2);

Yestim  = Ygrid(idx, :);
Xestim  = Xgrid(idx, :);