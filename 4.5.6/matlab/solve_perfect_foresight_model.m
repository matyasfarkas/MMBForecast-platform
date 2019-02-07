function [flag,endo_simul,err] = solve_perfect_foresight_model(endo_simul,exo_simul,pfm)

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

flag = 0;
err = 0;
stop = 0;
nan_flag = 0;

model_dynamic = pfm.dynamic_model;

Y = endo_simul(:);

if pfm.verbose
    disp (['-----------------------------------------------------']) ;
    disp (['MODEL SIMULATION :']) ;
    fprintf('\n') ;
end

if pfm.use_bytecode
    [flag, endo_simul]=bytecode(Y, exo_simul, pfm.params);
    return
end

z = Y(find(pfm.lead_lag_incidence'));
[d1,jacobian] = model_dynamic(z,exo_simul,pfm.params,pfm.steady_state,2);

% Initialization of the jacobian of the stacked model.
A = sparse([],[],[],pfm.periods*pfm.ny,pfm.periods*pfm.ny,pfm.periods*nnz(jacobian));

% Initialization of the Newton residuals.
res = zeros(pfm.periods*pfm.ny,1);

h1 = clock;

% Newton loop.
for iter = 1:pfm.maxit_
    h2 = clock;
    i_rows = 1:pfm.ny;
    i_cols = find(pfm.lead_lag_incidence');
    i_cols_A = i_cols;
    % Fill the jacobian of the stacked model.
    for it = 2:(pfm.periods+1)
        [d1,jacobian] = model_dynamic(Y(i_cols),exo_simul,pfm.params,pfm.steady_state,it);
        if it == 2
            A(i_rows,pfm.i_cols_A1) = jacobian(:,pfm.i_cols_1);
        elseif it == pfm.periods+1
            A(i_rows,i_cols_A(pfm.i_cols_T)) = jacobian(:,pfm.i_cols_T);
        else
            A(i_rows,i_cols_A) = jacobian(:,pfm.i_cols_j);
        end
        res(i_rows) = d1;
        i_rows = i_rows + pfm.ny;
        i_cols = i_cols + pfm.ny;
        if it > 2
            i_cols_A = i_cols_A + pfm.ny;
        end
    end
    % Stop if Newton residuals are zero.
    err = max(abs(res));
    if err < pfm.tolerance
        stop = 1 ;
        if pfm.verbose
            fprintf('\n') ;
            disp([' Total time of simulation        :' num2str(etime(clock,h1))]) ;
            fprintf('\n') ;
            disp([' Convergency obtained.']) ;
            fprintf('\n') ;
        end
        flag = 0;% Convergency obtained.
        endo_simul = reshape(Y,pfm.ny,pfm.periods+2);
        break
    end
    % Compute the Newton step.
    dy = -A\res;
    if any(isnan(dy))
        nan_flag = 1;
        break
    end
    % Update the endogenous variables paths.
    Y(pfm.i_upd) =   Y(pfm.i_upd) + dy;
end

if ~stop
    if pfm.verbose
        fprintf('\n') ;
        disp(['     Total time of simulation        :' num2str(etime(clock,h1))]) ;
        fprintf('\n') ;
        disp(['WARNING : maximum number of iterations is reached (modify options_.simul.maxit).']) ;
        fprintf('\n') ;
    end
    flag = 1;% more iterations are needed.
    endo_simul = 1;
end
if nan_flag
    if pfm.verbose
        fprintf('\n') ;
        disp(['     Total time of simulation        :' num2str(etime(clock,h1))]) ;
        fprintf('\n') ;
        disp(['WARNING : NaNs!']) ;
        fprintf('\n') ;
    end
    flag = 1;
    endo_simul = 1;
end
if pfm.verbose
    disp (['-----------------------------------------------------']) ;
end