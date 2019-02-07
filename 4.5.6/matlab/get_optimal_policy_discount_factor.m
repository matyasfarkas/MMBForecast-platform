function discount_factor=get_optimal_policy_discount_factor(params,param_names)

%function discount_factor=get_optimal_policy_discount_factor(M)
%  get the value of Ramsey policy discount factor
%
% INPUTS
%   params:             (vector) value of parameters
%   param_names:        (char array) list of parameter names
%
% OUTPUTS
%   discount_factor     (double) discount factor
%
% SPECIAL REQUIREMENTS
%   none

% Copyright (C) 2007-2017 Dynare Team
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

discount_factor = params(find(strcmp('optimal_policy_discount_factor',cellstr(param_names))));
