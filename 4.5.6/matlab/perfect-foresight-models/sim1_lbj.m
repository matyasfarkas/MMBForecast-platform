function [endogenousvariables, info] = sim1_lbj(endogenousvariables, exogenousvariables, steadystate, M, options)

% Performs deterministic simulations with lead or lag on one period using the historical LBJ algorithm
%
% INPUTS
%   ...
% OUTPUTS
%   ...
% ALGORITHM
%   Laffargue, Boucekkine, Juillard (LBJ)
%   see Juillard (1996) Dynare: A program for the resolution and
%   simulation of dynamic models with forward variables through the use
%   of a relaxation algorithm. CEPREMAP. Couverture Orange. 9602.
%
% SPECIAL REQUIREMENTS
%   None.

% Copyright (C) 1996-2017 Dynare Team
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

lead_lag_incidence = M.lead_lag_incidence;

ny = size(endogenousvariables,1);
nyp = nnz(lead_lag_incidence(1,:));
nyf = nnz(lead_lag_incidence(3,:));
nrs = ny+nyp+nyf+1;
nrc = nyf+1;
iyf = find(lead_lag_incidence(3,:)>0);
iyp = find(lead_lag_incidence(1,:)>0);
isp = [1:nyp];
is = [nyp+1:ny+nyp];
isf = iyf+nyp;
isf1 = [nyp+ny+1:nyf+nyp+ny+1];
stop = false;
iz = [1:ny+nyp+nyf];

dynamicmodel = sprintf('%s_dynamic', M.fname);

verbose = options.verbosity;

if verbose
    printline(56)
    disp(['MODEL SIMULATION :'])
    skipline()
end

it_init = M.maximum_lag+1;
h1 = clock;

for iter = 1:options.simul.maxit
    h2 = clock;
    if ~options.terminal_condition
        c = zeros(ny*options.periods, nrc);
    else
        c = zeros(ny*(options.periods+1), nrc);
    end
    it_ = it_init;
    z = [endogenousvariables(iyp,it_-1) ; endogenousvariables(:,it_) ; endogenousvariables(iyf,it_+1)];
    [d1, jacobian] = feval(dynamicmodel, z, exogenousvariables, M.params, steadystate, it_);
    jacobian = [jacobian(:,iz), -d1];
    ic = [1:ny];
    icp = iyp;
    c (ic,:) = jacobian(:,is)\jacobian(:,isf1);
    for it_ = it_init+(1:options.periods-1)
        z = [endogenousvariables(iyp,it_-1) ; endogenousvariables(:,it_) ; endogenousvariables(iyf,it_+1)];
        [d1, jacobian] = feval(dynamicmodel, z, exogenousvariables, M.params, steadystate, it_);
        jacobian = [jacobian(:,iz), -d1];
        jacobian(:,[isf nrs]) = jacobian(:,[isf nrs])-jacobian(:,isp)*c(icp,:);
        ic = ic + ny;
        icp = icp + ny;
        c (ic,:) = jacobian(:,is)\jacobian(:,isf1);
    end
    if options.terminal_condition == 1
        s = eye(ny);
        s(:,isf) = s(:,isf)+c(ic,1:nyf);
        ic = ic + ny;
        c(ic,nrc) = s\c(ic,nrc);
        c = bksup1(c, ny, nrc, iyf, options.periods);
        c = reshape(c, ny, options.periods+1);
        endogenousvariables(:,it_init+(0:options.periods)) = endogenousvariables(:,it_init+(0:options.periods))+options.slowc*c;
    else
        c = bksup1(c, ny, nrc, iyf, options.periods);
        c = reshape(c, ny, options.periods);
        endogenousvariables(:,it_init+(0:options.periods-1)) = endogenousvariables(:,it_init+(0:options.periods-1))+options.slowc*c;
    end
    err = max(max(abs(c./options.scalv')));
    if verbose
        str = sprintf('Iter: %s,\t err. = %s, \t time = %s', num2str(iter), num2str(err), num2str(etime(clock, h2)));
        disp(str);
    end
    if err < options.dynatol.f
        stop = true;
        if verbose
            skipline()
            disp(sprintf('Total time of simulation: %s', num2str(etime(clock,h1))))
        end
        info.status = 1;% Convergency obtained.
        info.error = err;
        info.iterations = iter;
        break
    end
end

if ~stop
    if verbose
        disp(sprintf('Total time of simulation: %s.', num2str(etime(clock,h1))))
        disp('Maximum number of iterations is reached (modify option maxit).')
    end
    info.status = 0;% more iterations are needed.
    info.error = err;
    info.errors = c/abs(err);
    info.iterations = options.simul.maxit;
end

if verbose
    if stop
        printline(56)
    else
        printline(62)
    end
    skipline()
end
