function [x,info] = dynare_solve_block_or_bytecode(y, exo, params, options, M)
% Copyright (C) 2010-2017 Dynare Team
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

info = 0;
x = y;
if options.block && ~options.bytecode
    for b = 1:length(M.block_structure_stat.block)
        ss = x;
        if M.block_structure_stat.block(b).Simulation_Type ~= 1 && ...
                M.block_structure_stat.block(b).Simulation_Type ~= 2
            if options.solve_algo <= 4
                [y, check] = dynare_solve('block_mfs_steadystate', ...
                                          ss(M.block_structure_stat.block(b).variable), ...
                                          options, b, ss, exo, params, M);
                if check ~= 0
                    %                    error(['STEADY: convergence
                    %                    problems in block ' int2str(b)])
                    info = 1;
                    return
                end
                ss(M.block_structure_stat.block(b).variable) = y;
            else
                n = length(M.block_structure_stat.block(b).variable);
                [ss, check] = solve_one_boundary([M.fname '_static_' int2str(b)], ss, exo, ...
                                                 params, [], M.block_structure_stat.block(b).variable, n, 1, 0, b, 0, options.simul.maxit, ...
                                                 options.solve_tolf, ...
                                                 options.slowc, 0, options.solve_algo, 1, 0, 0,M,options);
                if check
                    info = 1;
                    return
                end
            end
        end
        [r, g1, x] = feval([M.fname '_static'], b, ss, ...
                           exo, params);
    end
elseif options.bytecode
    if options.solve_algo > 4
        [check, x] = bytecode('static', x, exo, params);
        %        mexErrCheck('bytecode', check);
        if check
            info = 1;
            return
        end
    elseif options.block
        for b = 1:length(M.block_structure_stat.block)
            if M.block_structure_stat.block(b).Simulation_Type ~= 1 && ...
                    M.block_structure_stat.block(b).Simulation_Type ~= 2
                [y, check] = dynare_solve('block_bytecode_mfs_steadystate', ...
                                          x(M.block_structure_stat ...
                                            .block(b).variable), ...
                                          options, b, x, exo, params, M);
                if check
                    %                    error(['STEADY: convergence problems in block '
                    %                    int2str(b)])
                    info = 1;
                    return
                end
                x(M.block_structure_stat.block(b).variable) = y;
            else
                [chk, nulldev, nulldev1, x] = bytecode( x, exo, params, ...
                                                        x, 1, x, 'evaluate', 'static', ...
                                                        ['block = ' int2str(b)]);
                if chk
                    info = 1;
                    return
                end
            end
        end
    else
        [x, check] = dynare_solve('bytecode_steadystate', y, ...
                                  options, exo, params);
        if check
            %            error('STEADY: convergence problems')
            info = 1;
            return
        end
    end
end
