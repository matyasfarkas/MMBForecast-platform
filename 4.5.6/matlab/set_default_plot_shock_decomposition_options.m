function options = set_default_plot_shock_decomposition_options(options)
%function options = set_default_plot_shock_decomposition_options(options)
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

options.plot_shock_decomp.use_shock_groups = '';
options.plot_shock_decomp.colormap = '';
options.plot_shock_decomp.nodisplay = 0;
options.plot_shock_decomp.graph_format = 'eps';
options.plot_shock_decomp.detail_plot = 0;
options.plot_shock_decomp.interactive = 0;
options.plot_shock_decomp.screen_shocks = 0;
options.plot_shock_decomp.steadystate = 0;
options.plot_shock_decomp.type = '';
options.plot_shock_decomp.fig_name = '';
options.plot_shock_decomp.write_xls = 0;
options.plot_shock_decomp.realtime = 0; % 0 is standard; 1 is realtime
                                        % (pool/vintage); 2 is conditional
                                        % (pool/vintage); 3 is forecast
                                        % (pool/vintage)
options.plot_shock_decomp.vintage = 0; % 0 pool realtime/conditional; int:
                                       % forecast/conditional shock
                                       % decompositions
options.plot_shock_decomp.plot_init_date = [];
options.plot_shock_decomp.plot_end_date = [];
end
