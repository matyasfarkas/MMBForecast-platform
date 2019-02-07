function check_input_arguments(DynareOptions, DynareModel, DynareResults)
%function check_input_arguments(DynareOptions, DynareModel, DynareResults)
%Conducts checks for inconsistent/missing inputs to deterministic
%simulations

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

if DynareOptions.stack_solve_algo < 0 || DynareOptions.stack_solve_algo > 7
    error('perfect_foresight_solver:ArgCheck','PERFECT_FORESIGHT_SOLVER: stack_solve_algo must be between 0 and 7')
end

if ~DynareOptions.block && ~DynareOptions.bytecode && DynareOptions.stack_solve_algo ~= 0 ...
        && DynareOptions.stack_solve_algo ~= 6 && DynareOptions.stack_solve_algo ~= 7
    error('perfect_foresight_solver:ArgCheck','PERFECT_FORESIGHT_SOLVER: you must use stack_solve_algo=0 or stack_solve_algo=6 when not using block nor bytecode option')
end

if DynareOptions.block && ~DynareOptions.bytecode && DynareOptions.stack_solve_algo == 5
    error('perfect_foresight_solver:ArgCheck','PERFECT_FORESIGHT_SOLVER: you can''t use stack_solve_algo = 5 without bytecode option')
end

if (DynareOptions.block || DynareOptions.bytecode) && DynareOptions.stack_solve_algo == 6
    error('perfect_foresight_solver:ArgCheck','PERFECT_FORESIGHT_SOLVER: you can''t use stack_solve_algo = 6 with block or bytecode option')
end

if isoctave && DynareOptions.stack_solve_algo == 2
    error('perfect_foresight_solver:ArgCheck','PERFECT_FORESIGHT_SOLVER: you can''t use stack_solve_algo = 2 under Octave')
end


if isempty(DynareResults.endo_simul) || any(size(DynareResults.endo_simul) ~= [ DynareModel.endo_nbr, DynareModel.maximum_lag+DynareOptions.periods+DynareModel.maximum_lead ])

    if DynareOptions.initval_file
        fprintf('PERFECT_FORESIGHT_SOLVER: ''oo_.endo_simul'' has wrong size. Check whether your initval-file provides %d periods.',DynareModel.maximum_endo_lag+DynareOptions.periods+DynareModel.maximum_endo_lead)
        error('perfect_foresight_solver:ArgCheck','PERFECT_FORESIGHT_SOLVER: ''oo_.endo_simul'' has wrong size. Did you run ''perfect_foresight_setup'' ?')
    else
        error('perfect_foresight_solver:ArgCheck','PERFECT_FORESIGHT_SOLVER: ''oo_.endo_simul'' has wrong size. Did you run ''perfect_foresight_setup'' ?')
    end
end

if (DynareModel.exo_nbr > 0) && ...
        (isempty(DynareResults.exo_simul) || any(size(DynareResults.exo_simul) ~= [ DynareModel.maximum_lag+DynareOptions.periods+DynareModel.maximum_lead, DynareModel.exo_nbr ]))
    if DynareOptions.initval_file
        fprintf('PERFECT_FORESIGHT_SOLVER: ''oo_.exo_simul'' has wrong size. Check whether your initval-file provides %d periods.',DynareModel.maximum_endo_lag+DynareOptions.periods+DynareModel.maximum_endo_lead)
        error('perfect_foresight_solver:ArgCheck','PERFECT_FORESIGHT_SOLVER: ''oo_.exo_simul'' has wrong size.')
    else
        error('perfect_foresight_solver:ArgCheck','PERFECT_FORESIGHT_SOLVER: ''oo_.exo_simul'' has wrong size. Did you run ''perfect_foresight_setup'' ?')
    end
end
