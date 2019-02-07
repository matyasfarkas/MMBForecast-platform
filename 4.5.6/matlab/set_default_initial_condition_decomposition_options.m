function options = set_default_initial_condition_decomposition_options(options)
%function options = set_default_initial_condition_decomposition_options(options)
% sets the default options for prior_shock_decomposition
%
% INPUTS
%    options
%
% OUTPUTS
%    options
%
% SPECIAL REQUIREMENTS
%    none

% Copyright (C) 2017 Dynare Team
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

options.initial_condition_decomp.detail_plot = 0;
options.initial_condition_decomp.steadystate = 0;
options.initial_condition_decomp.write_xls = 0;
options.initial_condition_decomp.type = '';
options.initial_condition_decomp.plot_init_date = [];
options.initial_condition_decomp.plot_end_date = [];
end
