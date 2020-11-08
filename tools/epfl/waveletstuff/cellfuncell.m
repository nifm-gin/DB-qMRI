function c = cellfuncell(fun, c)
% c = cellfuncell(fun, c)
%
% Replacement for Matlab's cellfun routine that
% does not require a function that returns an integer.
%
% (c) Cedric Vonesch, 20077.10.19

for n = 1:numel(c)
	c{n} = feval(fun, c{n});
end