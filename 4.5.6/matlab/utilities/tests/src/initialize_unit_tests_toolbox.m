% Copyright (C) 2014-2017 Dynare Team
%
% This file is part of Dynare (m-unit-tests module).
%
% Dynare is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% Dynare's m-unit-tests module is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
% or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
% more details.
%
% You should have received a copy of the GNU General Public License
% along with Dynare.  If not, see <http://www.gnu.org/licenses/>.

% Get the path to the m-unit-tests/src folder.
unit_tests_src_root = strrep(which('initialize_unit_tests_toolbox'),'initialize_unit_tests_toolbox.m','');

if ~exist('isoctave','file')
    addpath([unit_tests_src_root '/missing/isoctave'])
end

if ~exist('skipline','file')
    addpath([unit_tests_src_root '/missing/skipline'])
end
