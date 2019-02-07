function planner_objective_value = evaluate_planner_objective(M,options,oo)

%function oo1 = evaluate_planner_objective(dr,M,oo,options)
%  computes value of planner objective function
%
% INPUTS
%   M:        (structure) model description
%   options:  (structure) options
%   oo:       (structure) output results
%
% SPECIAL REQUIREMENTS
%   none

% Copyright (C) 2007-2017 Dynare Team
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

dr = oo.dr;
exo_nbr = M.exo_nbr;
nstatic = M.nstatic;
nspred = M.nspred;
if nspred > 180
    disp(' ')
    disp(['WARNING in evaluate_planner_objective: model too large, can''t evaluate planner ' ...
          'objective'])
    planner_objective_value = NaN;
    return
end
beta = get_optimal_policy_discount_factor(M.params,M.param_names);

Gy = dr.ghx(nstatic+(1:nspred),:);
Gu = dr.ghu(nstatic+(1:nspred),:);
gy(dr.order_var,:) = dr.ghx;
gu(dr.order_var,:) = dr.ghu;


ys = oo.dr.ys;

[U,Uy,Uyy] = feval([M.fname '_objective_static'],ys,zeros(1,exo_nbr), ...
                   M.params);
%second order terms
Uyy = full(Uyy);

[Uyygygy, err] = A_times_B_kronecker_C(Uyy,gy,gy,options.threads.kronecker.A_times_B_kronecker_C);
mexErrCheck('A_times_B_kronecker_C', err);
[Uyygugu, err] = A_times_B_kronecker_C(Uyy,gu,gu,options.threads.kronecker.A_times_B_kronecker_C);
mexErrCheck('A_times_B_kronecker_C', err);
[Uyygygu, err] = A_times_B_kronecker_C(Uyy,gy,gu,options.threads.kronecker.A_times_B_kronecker_C);
mexErrCheck('A_times_B_kronecker_C', err);

Wbar =U/(1-beta); %steady state welfare
Wy = Uy*gy/(eye(nspred)-beta*Gy);
Wu = Uy*gu+beta*Wy*Gu;
Wyy = Uyygygy/(eye(nspred*nspred)-beta*kron(Gy,Gy));
[Wyygugu, err] = A_times_B_kronecker_C(Wyy,Gu,Gu,options.threads.kronecker.A_times_B_kronecker_C);
mexErrCheck('A_times_B_kronecker_C', err);
[Wyygygu,err] = A_times_B_kronecker_C(Wyy,Gy,Gu,options.threads.kronecker.A_times_B_kronecker_C);
mexErrCheck('A_times_B_kronecker_C', err);
Wuu = Uyygugu+beta*Wyygugu;
Wyu = Uyygygu+beta*Wyygygu;
Wss = beta*Wuu*M.Sigma_e(:)/(1-beta); % at period 0, we are in steady state, so the deviation term only starts in period 1, thus the beta in front

% initialize yhat1 at the steady state
yhat1 = oo.steady_state;
if options.ramsey_policy
    % initialize le Lagrange multipliers to 0 in yhat2
    yhat2 = zeros(M.endo_nbr,1);
    yhat2(1:M.orig_endo_nbr) = oo.steady_state(1:M.orig_endo_nbr);
end
if ~isempty(M.endo_histval)
    % initialize endogenous state variable to histval if necessary
    yhat1(1:M.orig_endo_nbr) = M.endo_histval(1:M.orig_endo_nbr);
    if options.ramsey_policy
        yhat2(1:M.orig_endo_nbr) = M.endo_histval(1:M.orig_endo_nbr);
    end
end
yhat1 = yhat1(dr.order_var(nstatic+(1:nspred)),1)-dr.ys(dr.order_var(nstatic+(1:nspred)));
u = oo.exo_simul(1,:)';

[Wyyyhatyhat1, err] = A_times_B_kronecker_C(Wyy,yhat1,yhat1,options.threads.kronecker.A_times_B_kronecker_C);
mexErrCheck('A_times_B_kronecker_C', err);
[Wuuuu, err] = A_times_B_kronecker_C(Wuu,u,u,options.threads.kronecker.A_times_B_kronecker_C);
mexErrCheck('A_times_B_kronecker_C', err);
[Wyuyhatu1, err] = A_times_B_kronecker_C(Wyu,yhat1,u,options.threads.kronecker.A_times_B_kronecker_C);
mexErrCheck('A_times_B_kronecker_C', err);
planner_objective_value(1) = Wbar+Wy*yhat1+Wu*u+Wyuyhatu1 ...
    + 0.5*(Wyyyhatyhat1 + Wuuuu+Wss);
if options.ramsey_policy
    yhat2 = yhat2(dr.order_var(nstatic+(1:nspred)),1)-dr.ys(dr.order_var(nstatic+(1:nspred)));
    [Wyyyhatyhat2, err] = A_times_B_kronecker_C(Wyy,yhat2,yhat2,options.threads.kronecker.A_times_B_kronecker_C);
    mexErrCheck('A_times_B_kronecker_C', err);
    [Wyuyhatu2, err] = A_times_B_kronecker_C(Wyu,yhat2,u,options.threads.kronecker.A_times_B_kronecker_C);
    mexErrCheck('A_times_B_kronecker_C', err);
    planner_objective_value(2) = Wbar+Wy*yhat2+Wu*u+Wyuyhatu2 ...
        + 0.5*(Wyyyhatyhat2 + Wuuuu+Wss);
end

if ~options.noprint
    skipline()
    disp('Approximated value of planner objective function')
    if options.ramsey_policy
        disp(['    - with initial Lagrange multipliers set to 0: ' ...
              num2str(planner_objective_value(2)) ])
        disp(['    - with initial Lagrange multipliers set to steady state: ' ...
              num2str(planner_objective_value(1)) ])
    elseif options.discretionary_policy
        fprintf('with discretionary policy: %10.8f',planner_objective_value(1))
    end
    skipline()
end
