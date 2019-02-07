function [residuals,check1,jacob] = evaluate_static_model(ys,exo_ss,params,M,options)

% function [residuals,check1,jacob] = evaluate_static_model(ys,exo_ss,params,M,options)
% Evaluates the static model
%
% INPUTS
%   ys                        vector           initial values used to compute the steady
%                                                 state
%   exo_ss                    vector           exogenous steady state
%   params                    vector           parameters
%   M                         struct           model structure
%   options                   struct           options
%
% OUTPUTS
%   residuals                 vector           residuals when ys is not
%                                              the steady state
%   check1                    scalar           error flag
%   jacob                     matrix           Jacobian of static model
%
% SPECIAL REQUIREMENTS
%   none

% Copyright (C) 2001-2017 Dynare Team
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

check1 = 0;
if options.bytecode
    [check1, residuals] = bytecode('evaluate','static',ys,...
                                   exo_ss, params, ys, 1);
    mexErrCheck('bytecode', check1);
else
    fh_static = str2func([M.fname '_static']);
    if options.block
        residuals = zeros(M.endo_nbr,1);
        for b = 1:length(M.block_structure_stat.block)
            mfsb = M.block_structure_stat.block(b).variable;
            % blocks that can be directly evaluated (mfsb is empty)
            % have zero residuals by construction
            if M.block_structure_stat.block(b).Simulation_Type ~= 1 && ...
                    M.block_structure_stat.block(b).Simulation_Type ~= 2
                residuals(mfsb) = feval(fh_static,b,ys,exo_ss,params);
            else
                %need to evaluate the recursive blocks to compute the
                %temporary terms
                feval(fh_static,b,ys,exo_ss,params);
            end
        end
    else
        residuals = feval(fh_static,ys,exo_ss,params);
    end
end
