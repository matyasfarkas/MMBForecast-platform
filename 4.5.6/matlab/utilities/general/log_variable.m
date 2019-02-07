function y = log_variable(ivar,x,M)
% returns the log of an endogenous variables except if
% it is a lagged/leaded exogenous variable
%
% INPUT
% ivar: vector of variable indices (in declaration order)
% x:    vector of variables value
% M:    model structure
%
% OUTPUT
% y:    log of the selected variables if there are not auxiliary variables
%       for lagged/leaded exogenous variables
%

% Copyright (C) 2011-2017 Dynare Team
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

orig_endo_nbr = M.orig_endo_nbr;
aux_vars = M.aux_vars;

y = zeros(length(ivar),1);
for i = 1:length(ivar)
    % Does ivar(i) refer to a lag/lead exogenous variable?
    if ivar(i) > orig_endo_nbr
        av = aux_vars(ivar(i) - orig_endo_nbr);
        if av.type == 2 || av.type == 3
            if av.endo_index ~= ivar(i)
                error(['This case shouldn''t happen. Please, contact Dynare ' ...
                       'developers'])
            else
                y(i) =  x(ivar(i));
            end
        else
            y(i) = log(x(ivar(i)));
        end
    else
        y(i) = log(x(ivar(i)));
    end
end