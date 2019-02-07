function err = evaluate_max_dynamic_residual(model_dynamic, Y, exogenous_variables, params, steady_state, periods, ny, max_lag, lead_lag_incidence)

% Copyright (C) 2013-2017 Dynare Team
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

i_rows = 1:ny;
i_cols = find(lead_lag_incidence');

err = 0;

for it = (max_lag+1):(max_lag+periods)
    d = model_dynamic(Y(i_cols), exogenous_variables, params, steady_state, it);
    i_rows = i_rows + ny;
    i_cols = i_cols + ny;
    r = max(abs(d));
    if r>err
        err = r;
    end
end