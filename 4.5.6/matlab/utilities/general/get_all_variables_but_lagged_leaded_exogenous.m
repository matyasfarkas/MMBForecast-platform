function k = get_all_variables_but_lagged_leaded_exogenous(M)
% returns indices of all endogenous variables in declaration order except
% lagged and leaded exogenous
%
% INPUT
% M: model structure
%
% OUTPUT
% k: vector of variable indices
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

if isempty(M.aux_vars)
    k = 1:M.endo_nbr;
else
    type = [M.aux_vars.type];
    k = [1:M.orig_endo_nbr, M.orig_endo_nbr ...
         + find((type ~= 2) & (type ~= 3))];
end