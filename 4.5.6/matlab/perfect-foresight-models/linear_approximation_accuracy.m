function err = linear_approximation_accuracy(options_, M_, oo_)
% Evaluates the accuracy of the linear approximation when solving perfect foresight models, by
% reporting the max absolute value of the dynamic residuals.
%
% INPUTS
% - options_ [struct] contains various options.
% - M_       [struct] contains a description of the model.
% - oo_      [struct] contains results.
%
% OUTPUTS
% - err      [double] n*1 vector, evaluation of the accuracy (n is the number of equations).

% Copyright (C) 2015-2017 Dynare Team
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

lead_lag_incidence = M_.lead_lag_incidence;

ny = M_.endo_nbr;

maximum_lag = M_.maximum_lag;

periods = options_.periods;
steady_state = oo_.steady_state;
params = M_.params;
endo_simul = oo_.endo_simul;
exo_simul = oo_.exo_simul;

model_dynamic = str2func([M_.fname,'_dynamic']);

residuals = zeros(ny,periods);

Y = endo_simul(:);

i_cols = find(lead_lag_incidence')+(maximum_lag-1)*ny;

for it = (maximum_lag+1):(maximum_lag+periods)
    residuals(:,it-1) = model_dynamic(Y(i_cols), exo_simul, params, steady_state,it);
    i_cols = i_cols + ny;
end

err = transpose(max(abs(transpose(residuals))));