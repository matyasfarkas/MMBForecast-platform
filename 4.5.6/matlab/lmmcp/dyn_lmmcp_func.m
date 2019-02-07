function [F,A] = dyn_lmmcp_func(x, model_dynamic, Y0, YT, exo_simul, params, ...
                                steady_state, periods, ny, lead_lag_incidence, ...
                                i_cols_A1, i_cols_1, i_cols_T, i_cols_j, ...
                                nnzA,eq_index)

% Copyright (C) 2014-2017 Dynare Team
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

Y = [Y0; x; YT];

F = zeros(periods*ny,1);
if nargout == 2
    A = sparse([],[],[],periods*ny,periods*ny,periods*nnzA);
end

i_rows = 1:ny;
i_cols = find(lead_lag_incidence');
i_cols_A = i_cols;

for it = 2:(periods+1)

    [res,jacobian] = model_dynamic(Y(i_cols),exo_simul, params, ...
                                   steady_state,it);
    F(i_rows) = res(eq_index);

    if nargout == 2
        if it == 2
            A(i_rows,i_cols_A1) = jacobian(eq_index,i_cols_1);
        elseif it == periods+1
            A(i_rows,i_cols_A(i_cols_T)) = jacobian(eq_index,i_cols_T);
        else
            A(i_rows,i_cols_A) = jacobian(eq_index,i_cols_j);
        end
    end

    i_rows = i_rows + ny;
    i_cols = i_cols + ny;
    if nargout == 2 && it > 2
        i_cols_A = i_cols_A + ny;
    end
end
