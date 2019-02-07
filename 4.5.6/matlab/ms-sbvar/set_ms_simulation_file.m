function options_=set_ms_simulation_file(options_)
% function options_=set_ms_simulation_file(options_)
% Set options_.ms.mh_file based on user input
%
% INPUTS
%    options_:    (struct)    options
%
% OUTPUTS
%    options_:    (struct)    options
%
% SPECIAL REQUIREMENTS
%    none

% Copyright (C) 2011-2012 Dynare Team
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

if isempty(options_.ms.simulation_file_tag)
    options_.ms.simulation_file_tag = options_.ms.file_tag;
end
options_.ms.mh_file = ['simulation_' options_.ms.simulation_file_tag '.out'];

if ~exist(options_.ms.mh_file,'file')
    error(['ERROR: Could not find Metropolis Hastings file: ' options_.ms.mh_file]);
end
end
