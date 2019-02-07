function [info_convergence, endo_simul] = extended_path_homotopy(endo_simul, exo_simul, M, options, oo, pfm, ep, order, algo, method, debug)

% Copyright (C) 2016-2017 Dynare Team
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


endo_simul0 = endo_simul;
if ismember(method, [1, 2])
    noconvergence = true;
    iteration = 0;
    weight = .1;
    maxiter = 100;
    increase_flag = false;
    increase_factor = 1.2;
    decrease_factor = 1.1;
    state = false(5,1);
    oldweight = weight;
    while noconvergence
        iteration = iteration + 1;
        oo.endo_simul = endo_simul;
        oo.endo_simul(:,1) = oo.steady_state + weight*(endo_simul0(:,1) - oo.steady_state);
        oo.exo_simul = bsxfun(@plus, weight*exo_simul, (1-weight)*transpose(oo.exo_steady_state));
        if order==0
            tmp = perfect_foresight_solver_core(M, options, oo);
        else
            switch(algo)
              case 0
                [flag, tmp.endo_simul] = ...
                    solve_stochastic_perfect_foresight_model(endo_simul, exo_simul, pfm, ep.stochastic.quadrature.nodes, ep.stochastic.order);
              case 1
                [flag, tmp.endo_simul] = ...
                    solve_stochastic_perfect_foresight_model_1(endo_simul, exo_simul, options, pfm, ep.stochastic.order);
            end
        end
        if isequal(order, 0)
            % Logical variable flag is false iff the solver fails.
            flag = ~tmp.deterministic_simulation.status;
        else
            % Fix convention issue on the value of flag.
            flag = ~flag;
        end
        if debug
            disp(sprintf('%s\t %1.8f\t %s',int2str(iteration),weight,int2str(flag)))
        end
        state(2:end) = state(1:end-1);
        state(1) = flag;
        if flag
            if isequal(weight, 1)
                noconvergence = false;
                break
            end
            if all(state)
                increase_factor = 1+(increase_factor-1)*1.1;
                state = false(size(state));
            end
            oldweight = weight;
            weight = min(weight*increase_factor, 1);
            increase_flag = true;
            endo_simul = tmp.endo_simul;
        else
            if increase_flag
                weight = oldweight + (weight-oldweight)/100;
            else
                weight = min(weight/decrease_factor, 1);
            end
        end
        if iteration>maxiter
            break
        end
        if weight<1e-9
            break
        end
    end
    info_convergence = ~noconvergence;
end

if isequal(method, 3) || (isequal(method, 2) && noconvergence)
    if isequal(method, 2)
        endo_simul = endo_simul0;
    end
    weights = 0:(1/1000):1;
    noconvergence = true;
    index = 1;
    jndex = 0;
    nweights = length(weights);
    while noconvergence
        weight = weights(index);
        oo.endo_simul = endo_simul;
        oo.exo_simul = bsxfun(@plus, weight*exo_simul, (1-weight)*transpose(oo.exo_steady_state));
        if order==0
            tmp = perfect_foresight_solver_core(M, options, oo);
        else
            switch(algo)
              case 0
                [flag, tmp.endo_simul] = ...
                    solve_stochastic_perfect_foresight_model(endo_simul, exo_simul, pfm, ep.stochastic.quadrature.nodes, ep.stochastic.order);
              case 1
                [flag, tmp.endo_simul] = ...
                    solve_stochastic_perfect_foresight_model_1(endo_simul, exo_simul, options, pfm, ep.stochastic.order);
            end
        end
        if isequal(order, 0)
            % Logical variable flag is false iff the solver fails.
            flag = ~tmp.deterministic_simulation.status;
        else
            % Fix convention issue on the value of flag.
            flag = ~flag;
        end
        if debug
            disp(sprintf('%s\t %1.8f\t %s',int2str(index),weight,int2str(flag)))
        end
        if flag
            jndex = index;
            if isequal(weight, 1)
                noconvergence = false;
                continue
            end
            index = index+1;
            endo_simul = tmp.endo_simul;
        else
            break
        end
    end
    info_convergence = ~noconvergence;
end