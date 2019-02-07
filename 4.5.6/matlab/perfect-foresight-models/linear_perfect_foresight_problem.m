function [residuals,JJacobian] = linear_perfect_foresight_problem(y, dynamicjacobian, Y0, YT, ...
                                                  exo_simul, params, steady_state, ...
                                                  maximum_lag, T, ny, i_cols, ...
                                                  i_cols_J1, i_cols_1, i_cols_T, ...
                                                  i_cols_j,nnzJ,jendo,jexog)
% function [residuals,JJacobian] = perfect_foresight_problem(x, model_dynamic, Y0, YT,exo_simul,
% params, steady_state, maximum_lag, periods, ny, i_cols,i_cols_J1, i_cols_1,
% i_cols_T, i_cols_j, nnzA)
% computes the residuals and th Jacobian matrix
% for a perfect foresight problem over T periods.
%
% INPUTS
%   ...
% OUTPUTS
%   ...
% ALGORITHM
%   ...
%
% SPECIAL REQUIREMENTS
%   None.

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


YY = [Y0; y; YT];

residuals = zeros(T*ny,1);

z = zeros(columns(dynamicjacobian), 1);

if nargout == 2
    JJacobian = sparse([],[],[],T*ny,T*ny,T*nnzJ);
end

i_rows = 1:ny;
i_cols_J = i_cols;

for it = maximum_lag+(1:T)
    z(jendo) = YY(i_cols);
    z(jexog) = transpose(exo_simul(it,:));
    residuals(i_rows) = dynamicjacobian*z;
    if nargout == 2
        if it == maximum_lag+1
            JJacobian(i_rows,i_cols_J1) = dynamicjacobian(:,i_cols_1);
        elseif it == maximum_lag+T
            JJacobian(i_rows,i_cols_J(i_cols_T)) = dynamicjacobian(:,i_cols_T);
        else
            JJacobian(i_rows,i_cols_J) = dynamicjacobian(:,i_cols_j);
            i_cols_J = i_cols_J + ny;
        end
    end
    i_rows = i_rows + ny;
    i_cols = i_cols + ny;
end