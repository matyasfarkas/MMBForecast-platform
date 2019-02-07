function [M,oo,info,ip,ix,ixd]=homotopy3(values, step_nbr, M, options, oo)
% function homotopy3(values, step_nbr)
%
% Implements homotopy (mode 3) for steady-state computation.
% Tries first the most extreme values. If it fails to compute the steady
% state, the interval between initial and desired values is divided by two
% for each parameter. Every time that it is impossible to find a steady
% state, the previous interval is divided by two. When one succeed to find
% a steady state, the previous interval is multiplied by two.
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
%    step_nbr:      maximum number of steps to try before aborting
%    M              struct of model parameters
%    options        struct of options
%    oo             struct of outputs
%
% OUTPUTS
%    M              struct of model parameters
%    oo             struct of outputs
%    info           return status 0: OK, 1: failed
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

info = [];
tol = 1e-8;

nv = size(values,1);

ip = find(values(:,1) == 4); % Parameters
ix = find(values(:,1) == 1); % Exogenous
ixd = find(values(:,1) == 2); % Exogenous deterministic

if length([ip; ix; ixd]) ~= nv
    error('HOMOTOPY mode 3: incorrect variable types specified')
end

% Construct vector of starting values, using previously initialized values
% when initial value has not been given in homotopy_setup block
last_values = values(:,3);
ipn = find(values(:,1) == 4 & isnan(last_values));
last_values(ipn) = M.params(values(ipn, 2));
ixn = find(values(:,1) == 1 & isnan(last_values));
last_values(ixn) = oo.exo_steady_state(values(ixn, 2));
ixdn = find(values(:,1) == 2 & isnan(last_values));
last_values(ixdn) = oo.exo_det_steady_state(values(ixdn, 2));

targetvalues = values(:,4);

%if min(abs(targetvalues - last_values)) < tol
%    error('HOMOTOPY mode 3: distance between initial and final values should be at least %e for all variables', tol)
%end
iplus = find(targetvalues > last_values);
iminus = find(targetvalues < last_values);

curvalues = last_values;
inc = (targetvalues-last_values)/2;
kplus = [];
kminus = [];
last_values = [];

disp('HOMOTOPY mode 3: launching solver at initial point...')

iter = 1;
while iter <= step_nbr

    M.params(values(ip,2)) = curvalues(ip);
    oo.exo_steady_state(values(ix,2)) = curvalues(ix);
    oo.exo_det_steady_state(values(ixd,2)) = curvalues(ixd);

    old_ss = oo.steady_state;

    [steady_state,params,info] = steady_(M,options,oo);
    if info(1) == 0
        oo.steady_state = steady_state;
        M.params = params;
        if length([kplus; kminus]) == nv
            return
        end
        if iter == 1
            disp('HOMOTOPY mode 3: successful step, now jumping to final point...')
        else
            disp('HOMOTOPY mode 3: successful step, now multiplying increment by 2...')
        end
        last_values = curvalues;
        old_params = params;
        old_exo_steady_state = oo.exo_steady_state;
        old_exo_det_steady_state = oo.exo_det_steady_state;
        inc = 2*inc;
    elseif iter == 1
        error('HOMOTOPY mode 3: can''t solve the model at 1st iteration')
    else
        disp('HOMOTOPY mode 3: failed step, now dividing increment by 2...')
        inc = inc/2;
        oo.steady_state = old_ss;
    end

    curvalues = last_values + inc;
    kplus = find(curvalues(iplus) >= targetvalues(iplus));
    curvalues(iplus(kplus)) = targetvalues(iplus(kplus));
    kminus = find(curvalues(iminus) <= targetvalues(iminus));
    curvalues(iminus(kminus)) = targetvalues(iminus(kminus));

    if max(abs(inc)) < tol
        disp('HOMOTOPY mode 3: failed, increment has become too small')
        M.params = old_params;
        oo.exo_steady_state = old_exo_steady_state;
        oo.exo_det_steady_state = old_exo_det_steady_state;
        return
    end

    iter = iter + 1;
end
disp('HOMOTOPY mode 3: failed, maximum iterations reached')
M.params = old_params;
oo.exo_steady_state = old_exo_steady_state;
oo.exo_det_steady_state = old_exo_det_steady_state;
