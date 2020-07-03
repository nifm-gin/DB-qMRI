function [Yestim, score, Xestim] = EstimateParametersFromGrid(Xobs, Xgrid, Ygrid, verb)

narginchk(3, 4)
if nargin == 3, verb = 0; end

% normalization
Xobs_normalized = Xobs ./ vecnorm(Xobs,2,2);
Xgrid_normalized = Xgrid ./ vecnorm(Xgrid,2,2);

% dot-product/scalar product comparison
[~, score] = max(abs(Xobs_normalized * Xgrid_normalized'), [], 2);

Yestim  = Ygrid(score,:);
Xestim  = Xgrid(score,:);