function [endogenousvariables, info] = sim1_purely_backward(endogenousvariables, exogenousvariables, steadystate, M, options)

% Performs deterministic simulation of a purely backward model

% Copyright (C) 2012-2017 Dynare Team
%
% This file is part of Dynare.
%
% Dynare is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% Dynare is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with Dynare.  If not, see <http://www.gnu.org/licenses/>.

if size(M.lead_lag_incidence,1) > 1
    ny0 = nnz(M.lead_lag_incidence(2,:));    % Number of variables at current period
    nyb = nnz(M.lead_lag_incidence(1,:));    % Number of variables at previous period
    iyb = find(M.lead_lag_incidence(1,:)>0); % Indices of variables at previous period
else
    ny0 = nnz(M.lead_lag_incidence(1,:));    % Number of variables at current period
    nyb = 0;
    iyb = [];
end


if ny0 ~= M.endo_nbr
    error('All endogenous variables must appear at the current period!')
end

dynamicmodel = str2func([M.fname,'_dynamic']);

info.status = 1;

for it = 2:options.periods+1
    yb = endogenousvariables(:,it-1);        % Values at previous period, also used as guess value for current period
    yb1 = yb(iyb);
    [tmp, check] = solve1(dynamicmodel, [yb1; yb], 1:M.endo_nbr, nyb+1:nyb+M.endo_nbr, ...
                          1, options.gstep, options.solve_tolf, options.solve_tolx, ...
                          options.simul.maxit, options.debug, exogenousvariables, ...
                          M.params, steadystate, it+M.maximum_lag-1);
    if check
        info.status = 0;
    end
    endogenousvariables(:,it) = tmp(nyb+1:nyb+M.endo_nbr);
end