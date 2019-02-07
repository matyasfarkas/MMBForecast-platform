function clean_ms_simulation_files(file_tag)
% function clean_ms_simulation_files(file_tag)
% removes MS simulation files
%
% INPUTS
%    file_tag: string indicating tag to use when deleting files
%
% OUTPUTS
%    none
%
% SPECIAL REQUIREMENTS
%    none

% Copyright (C) 2011-2013 Dynare Team
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

delete_if_exists(['simulation_' file_tag '.out']);
delete_if_exists(['simulation_info_' file_tag '.out']);
delete_if_exists(['draws_test_' file_tag '.out']);
delete_if_exists(['draws_header_' file_tag '.out']);
end
