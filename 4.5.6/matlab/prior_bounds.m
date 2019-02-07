function bounds = prior_bounds(bayestopt, prior_trunc)

%@info:
%! @deftypefn {Function File} {@var{bounds} =} prior_bounds (@var{bayesopt},@var{option})
%! @anchor{prior_bounds}
%! @sp 1
%! Returns bounds for the prior densities. For each estimated parameter the lower and upper bounds
%! are such that the defined intervals contains a probability mass equal to 1-2*@var{option}.prior_trunc. The
%! default value for @var{option}.prior_trunc is 1e-10 (set in @ref{global_initialization}).
%! @sp 2
%! @strong{Inputs}
%! @sp 1
%! @table @ @var
%! @item bayestopt
%! Matlab's structure describing the prior distribution (initialized by @code{dynare}).
%! @item option
%! Matlab's structure describing the options (initialized by @code{dynare}).
%! @end table
%! @sp 2
%! @strong{Outputs}
%! @sp 1
%! @table @ @var
%! @item bounds
%! A structure with two fields lb and up (p*1 vectors of doubles, where p is the number of estimated parameters) for the lower and upper bounds.
%! @end table
%! @sp 2
%! @strong{This function is called by:}
%! @sp 1
%! @ref{get_prior_info}, @ref{dynare_estimation_1}, @ref{dynare_estimation_init}
%! @sp 2
%! @strong{This function calls:}
%! @sp 1
%! None.
%! @end deftypefn
%@eod:


% function bounds = prior_bounds(bayestopt)
% computes bounds for prior density.
%
% INPUTS
%    bayestopt  [structure]  characterizing priors (shape, mean, p1..p4)
%
% OUTPUTS
%    bounds     [double]      structure specifying prior bounds (lb and ub fields)
%
% SPECIAL REQUIREMENTS
%    none

% Copyright (C) 2003-2017 Dynare Team
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

pshape = bayestopt.pshape;
p3 = bayestopt.p3;
p4 = bayestopt.p4;
p6 = bayestopt.p6;
p7 = bayestopt.p7;

bounds.lb = zeros(length(p6),1);
bounds.ub = zeros(length(p6),1);

for i=1:length(p6)
    switch pshape(i)
      case 1
        if prior_trunc == 0
            bounds.lb(i) = p3(i);
            bounds.ub(i) = p4(i);
        else
            bounds.lb(i) = betainv(prior_trunc,p6(i),p7(i))*(p4(i)-p3(i))+p3(i);
            bounds.ub(i) = betainv(1-prior_trunc,p6(i),p7(i))*(p4(i)-p3(i))+p3(i);
        end
      case 2
        if prior_trunc == 0
            bounds.lb(i) = p3(i);
            bounds.ub(i) = Inf;
        else
            try
                bounds.lb(i) = gaminv(prior_trunc,p6(i),p7(i))+p3(i);
                bounds.ub(i) = gaminv(1-prior_trunc,p6(i),p7(i))+p3(i);
            catch
                % Workaround for ticket #161
                if isoctave
                    error(['Due to a bug in Octave, you must choose other values for mean and/or variance of your prior on ' bayestopt.name{i} ', or use another shape'])
                else
                    rethrow(lasterror)
                end
            end
        end
      case 3
        if prior_trunc == 0
            bounds.lb(i) = -Inf;
            bounds.ub(i) = Inf;
        else
            bounds.lb(i) = norminv(prior_trunc,p6(i),p7(i));
            bounds.ub(i) = norminv(1-prior_trunc,p6(i),p7(i));
        end
      case 4
        if prior_trunc == 0
            bounds.lb(i) = p3(i);
            bounds.ub(i) = Inf;
        else
            try
                bounds.lb(i) = 1/sqrt(gaminv(1-prior_trunc, p7(i)/2, 2/p6(i)))+p3(i);
                bounds.ub(i) = 1/sqrt(gaminv(prior_trunc, p7(i)/2, 2/p6(i)))+p3(i);
            catch
                % Workaround for ticket #161
                if isoctave
                    error(['Due to a bug in Octave, you must choose other values for mean and/or variance of your prior on ' bayestopt.name{i} ', or use another shape'])
                else
                    rethrow(lasterror)
                end
            end
        end
      case 5
        if prior_trunc == 0
            bounds.lb(i) = p6(i);
            bounds.ub(i) = p7(i);
        else
            bounds.lb(i) = p6(i)+(p7(i)-p6(i))*prior_trunc;
            bounds.ub(i) = p7(i)-(p7(i)-p6(i))*prior_trunc;
        end
      case 6
        if prior_trunc == 0
            bounds.lb(i) = p3(i);
            bounds.ub(i) = Inf;
        else
            try
                bounds.lb(i) = 1/gaminv(1-prior_trunc, p7(i)/2, 2/p6(i))+p3(i);
                bounds.ub(i) = 1/gaminv(prior_trunc, p7(i)/2, 2/p6(i))+ p3(i);
            catch
                % Workaround for ticket #161
                if isoctave
                    error(['Due to a bug in Octave, you must choose other values for mean and/or variance of your prior on ' bayestopt.name{i} ', or use another shape'])
                else
                    rethrow(lasterror)
                end
            end
        end
      case 8
        if prior_trunc == 0
            bounds.lb(i) = p3(i);
            bounds.ub(i) = Inf;
        else
            bounds.lb(i) = p3(i)+wblinv(prior_trunc,p6(i),p7(i));
            bounds.ub(i) = p3(i)+wblinv(1-prior_trunc,p6(i),p7(i));
        end
      otherwise
        error(sprintf('prior_bounds: unknown distribution shape (index %d, type %d)', i, pshape(i)));
    end
end