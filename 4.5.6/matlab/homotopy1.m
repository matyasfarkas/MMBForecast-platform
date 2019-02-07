function [M,oo,info,ip,ix,ixd] = homotopy1(values, step_nbr, M, options, oo)
% function homotopy1(values, step_nbr)
%
% Implements homotopy (mode 1) for steady-state computation.
% The multi-dimensional vector going from the set of initial values
% to the set of final values is divided in as many sub-vectors as
% there are steps, and the problem is solved as many times.
%
% INPUTS
%    values:        a matrix with 4 columns, representing the content of
%                   homotopy_setup block, with one variable per line.
%                   Column 1 is variable type (1 for exogenous, 2 for
%                   exogenous deterministic, 4 for parameters)
%                   Column 2 is symbol integer identifier.
%                   Column 3 is initial value, and column 4 is final value.
%                   Column 3 can contain NaNs, in which case previous
%                   initialization of variable will be used as initial value.
%    step_nbr:      number of steps for homotopy
%    M              struct of model parameters
%    options        struct of options
%    oo             struct of outputs
%
% OUTPUTS
%    M              struct of model parameters
%    oo             struct of outputs
%    ip             index of parameters
%    ix             index of exogenous variables
%    ixp            index of exogenous deterministic variables
%
% SPECIAL REQUIREMENTS
%    none

% Copyright (C) 2008-2017 Dynare Team
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

nv = size(values, 1);

ip = find(values(:,1) == 4); % Parameters
ix = find(values(:,1) == 1); % Exogenous
ixd = find(values(:,1) == 2); % Exogenous deterministic

if length([ip; ix; ixd]) ~= nv
    error('HOMOTOPY mode 1: incorrect variable types specified')
end

% Construct vector of starting values, using previously initialized values
% when initial value has not been given in homotopy_setup block
oldvalues = values(:,3);
ipn = find(values(:,1) == 4 & isnan(oldvalues));
oldvalues(ipn) = M.params(values(ipn, 2));
ixn = find(values(:,1) == 1 & isnan(oldvalues));
oldvalues(ixn) = oo.exo_steady_state(values(ixn, 2));
ixdn = find(values(:,1) == 2 & isnan(oldvalues));
oldvalues(ixdn) = oo.exo_det_steady_state(values(ixdn, 2));

points = zeros(nv, step_nbr+1);
for i = 1:nv
    if (oldvalues(i) ~= values(i, 4))
        points(i,:) = oldvalues(i):(values(i,4)-oldvalues(i))/step_nbr:values(i,4);
    else
        points(i,:) = values(i,4);
    end
end

for i=1:step_nbr+1
    disp([ 'HOMOTOPY mode 1: computing step ' int2str(i-1) '/' int2str(step_nbr) '...' ])
    old_params = M.params;
    old_exo = oo.exo_steady_state;
    old_exo_det = oo.exo_det_steady_state;
    M.params(values(ip,2)) = points(ip,i);
    oo.exo_steady_state(values(ix,2)) = points(ix,i);
    oo.exo_det_steady_state(values(ixd,2)) = points(ixd,i);

    [steady_state,M.params,info] = steady_(M,options,oo);
    if info(1) == 0
        % if homotopy step is not successful, current values of steady
        % state are not modified
        oo.steady_state = steady_state;
    else
        M.params = old_params;
        oo.exo_steady_state = old_exo;
        oo.exo_det_steady_state = old_exo_det;
        break
    end
end
